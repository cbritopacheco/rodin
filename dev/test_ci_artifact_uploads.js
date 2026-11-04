'use strict';

const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const os = require('node:os');
const path = require('node:path');
const {test} = require('node:test');
const {readInputs, transient, uploadArtifact} = require('../.github/actions/upload-artifacts/upload.js');
const clientRoot = process.env.ARTIFACT_CLIENT_ROOT || path.resolve('.github/actions/download-artifacts');
const glob = require(path.join(clientRoot, 'node_modules/@actions/glob'));
const digest = 'a'.repeat(64);
const network = () => Object.assign(new Error('Failed to CreateArtifact: ECONNRESET https://secret.example/?sig=private'), {code: 'ECONNRESET'});

async function fixture(t, overrides = {}) {
  const root = await fs.mkdtemp(path.join(os.tmpdir(), 'rodin-upload-'));
  t.after(() => fs.rm(root, {recursive: true, force: true}));
  await fs.mkdir(path.join(root, 'payload'));
  await fs.writeFile(path.join(root, 'payload', 'data.bin'), 'complete payload');
  const messages = [], outputs = {};
  const core = {info: () => {}, warning: text => messages.push(text), setOutput: (key, value) => {outputs[key] = value;}};
  const input = readInputs({UPLOAD_NAME: 'ci-build-test', UPLOAD_PATH: path.join(root, 'payload'),
    UPLOAD_OVERWRITE: 'true', UPLOAD_COMPRESSION: '0', UPLOAD_RETENTION: '7', UPLOAD_NO_FILES: 'error', ...overrides});
  return {root, messages, outputs, core, input};
}

test('CreateArtifact connection reset retries with the same immutable files and options', async t => {
  const f = await fixture(t);let calls = 0, deletions = 0;const waits = [];
  const client = {
    deleteArtifact: async name => {assert.equal(name, f.input.name);++deletions;},
    uploadArtifact: async (name, files, root, options) => {
      ++calls;
      assert.equal(name, f.input.name);
      assert.equal(root, path.join(f.root, 'payload'));
      assert.deepEqual(options, {compressionLevel: 0, retentionDays: 7});
      assert.deepEqual(files, [path.join(f.root, 'payload', 'data.bin')]);
      assert.equal(await fs.readFile(files[0], 'utf8'), 'complete payload');
      assert.deepEqual(f.outputs, {});
      if (calls === 1) throw network();
      return {id: 42, digest};
    }
  };
  await uploadArtifact(client, glob, f.input, f.core, async delay => waits.push(delay));
  assert.equal(calls, 2);assert.equal(deletions, 2);assert.deepEqual(waits, [2000]);
  assert.equal(f.outputs['artifact-id'], 42);assert.equal(f.outputs['artifact-digest'], digest);
  assert.equal(f.messages.length, 1);assert(!f.messages[0].includes('private'));
});

test('lost finalization response replaces a finalized artifact before retrying', async t => {
  const f = await fixture(t);let exists = false, calls = 0;
  const client = {
    deleteArtifact: async () => {exists = false;},
    uploadArtifact: async () => {
      assert.equal(exists, false);exists = true;
      if (++calls === 1) throw network();
      return {id: 43, digest};
    }
  };
  await uploadArtifact(client, glob, f.input, f.core, async () => {});
  assert.equal(calls, 2);assert.equal(f.outputs['artifact-id'], 43);
});

test('network exhaustion fails after exactly three attempts with no success outputs', async t => {
  const f = await fixture(t);let calls = 0;const waits = [];
  await assert.rejects(uploadArtifact({deleteArtifact: async () => {}, uploadArtifact: async () => {++calls;throw network();}},
    glob, f.input, f.core, async delay => waits.push(delay)), /ECONNRESET/);
  assert.equal(calls, 3);assert.deepEqual(waits, [2000, 4000]);assert.deepEqual(f.outputs, {});
});

test('permission, quota, missing files and cancellation do not retry', async t => {
  const f = await fixture(t);
  for (const error of [Object.assign(new Error('denied ECONNRESET'), {status: 403}),
    new Error('Artifact storage quota exceeded'), Object.assign(new Error('cancelled'), {name: 'AbortError'})]) {
    let calls = 0;
    await assert.rejects(uploadArtifact({deleteArtifact: async () => {}, uploadArtifact: async () => {++calls;throw error;}},
      glob, f.input, f.core, async () => assert.fail('Unexpected retry')));
    assert.equal(calls, 1);assert.deepEqual(f.outputs, {});
  }
  f.input.patterns = path.join(f.root, 'absent');
  await assert.rejects(uploadArtifact({uploadArtifact: async () => assert.fail('No files')}, glob, f.input, f.core), /no matching files/);
});

test('deletion retries transient errors and never ignores authorization errors', async t => {
  const f = await fixture(t);let calls = 0, uploads = 0;
  const client = {deleteArtifact: async () => {if (++calls === 1) throw network();},
    uploadArtifact: async () => {++uploads;return {id: 44, digest};}};
  await uploadArtifact(client, glob, f.input, f.core, async () => {});
  assert.equal(calls, 2);assert.equal(uploads, 1);
  client.deleteArtifact = async () => {throw Object.assign(new Error('forbidden'), {statusCode: 403});};
  await assert.rejects(uploadArtifact(client, glob, f.input, f.core, async () => assert.fail('Forbidden retry')), /forbidden/);
  assert.equal(uploads, 1);
});

test('exact missing-artifact deletion is allowed; non-overwrite never deletes', async t => {
  const f = await fixture(t);let deletes = 0;
  const client = {deleteArtifact: async () => {++deletes;throw Object.assign(new Error('absent'), {name: 'ArtifactNotFoundError'});},
    uploadArtifact: async () => ({id: 45, digest})};
  await uploadArtifact(client, glob, f.input, f.core);
  f.input.overwrite = false;
  await uploadArtifact(client, glob, f.input, f.core);
  assert.equal(deletes, 1);
});

test('glob exclusions, hidden files, multiple roots and single-file layout are preserved', async t => {
  const f = await fixture(t);const folder = path.join(f.root, 'payload');
  await fs.writeFile(path.join(folder, 'excluded.log'), 'exclude');
  await fs.writeFile(path.join(folder, '.hidden'), 'hidden');
  await fs.mkdir(path.join(f.root, 'second'));
  await fs.writeFile(path.join(f.root, 'second', 'other.txt'), 'other');
  const seen = [];
  const client = {deleteArtifact: async () => {}, uploadArtifact: async (name, files, root) => {
    seen.push({files: files.map(file => path.relative(f.root, file)).sort(), root});return {id: 46, digest};}};
  f.input.patterns = `${folder}\n!${folder}/excluded.log\n${f.root}/second`;
  await uploadArtifact(client, glob, f.input, f.core);
  assert.deepEqual(seen[0], {root: f.root, files: ['payload/data.bin', 'second/other.txt']});
  f.input.patterns = path.join(folder, 'data.bin');
  await uploadArtifact(client, glob, f.input, f.core);
  assert.equal(seen[1].root, folder);
  f.input.patterns = folder;f.input.includeHidden = true;
  await uploadArtifact(client, glob, f.input, f.core);
  assert(seen[2].files.includes('payload/.hidden'));
});

test('warn and ignore missing-file policies avoid any network calls', async t => {
  const f = await fixture(t);f.input.patterns = path.join(f.root, 'absent');
  for (const noFiles of ['warn', 'ignore']) {
    await uploadArtifact({uploadArtifact: async () => assert.fail('Optional absent artifact')}, glob,
      {...f.input, noFiles}, f.core);
  }
  assert.equal(f.messages.length, 1);assert.deepEqual(f.outputs, {});
});

test('only finalized uploads produce outputs, and option/status validation fails closed', async t => {
  const f = await fixture(t);
  await assert.rejects(uploadArtifact({deleteArtifact: async () => {}, uploadArtifact: async () => ({id: 0})},
    glob, f.input, f.core), /finalized/);
  assert.deepEqual(f.outputs, {});
  assert(transient(Object.assign(new Error('unavailable'), {statusCode: 503})));
  assert(transient(new Error('wrapped', {cause: network()})));
  assert(!transient(Object.assign(new Error('conflict'), {status: 409})));
  assert(transient(new Error('Failed to CreateArtifact: Failed request: (503) Service Unavailable')));
  assert(!transient(new Error('Failed to CreateArtifact: Received non-retryable error: Failed request: (403) Forbidden')));
  assert(transient(new Error('Upload progress stalled.')));
  for (const values of [{UPLOAD_COMPRESSION: '10'}, {UPLOAD_OVERWRITE: 'yes'}, {UPLOAD_RETENTION: '-1'}]) {
    assert.throws(() => readInputs({UPLOAD_NAME: 'x', UPLOAD_PATH: 'x', ...values}), /Invalid upload option/);
  }
});
