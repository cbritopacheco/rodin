'use strict';

// Real pinned-client RPC calls, ZIP streams and blob HTTP transfers against a
// local service; no GitHub account, repository token or remote artifact writes.
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const fs = require('node:fs/promises');
const http = require('node:http');
const os = require('node:os');
const path = require('node:path');
const {Readable} = require('node:stream');
const {test} = require('node:test');
const {uploadArtifact, readInputs} = require('../.github/actions/upload-artifacts/upload.js');
const root = process.env.ARTIFACT_CLIENT_ROOT || path.resolve('.github/actions/download-artifacts');
const {DefaultArtifactClient} = require(path.join(root, 'node_modules/@actions/artifact'));
const glob = require(path.join(root, 'node_modules/@actions/glob'));
const unzip = require(path.join(root, 'node_modules/unzip-stream'));

async function fixture(t, failure) {
  const folder = await fs.mkdtemp(path.join(os.tmpdir(), 'rodin-upload-stream-'));
  const payload = Buffer.alloc(65536, 29);
  await fs.writeFile(path.join(folder, 'payload.bin'), payload);
  const state = {creates: 0, finalizes: 0, deletes: 0, exists: false, blocks: new Map(), zip: null};
  let base;
  const server = http.createServer(async (request, response) => {
    const chunks = [];
    for await (const chunk of request) chunks.push(chunk);
    const body = Buffer.concat(chunks);
    const url = new URL(request.url, base);
    if (request.method === 'PUT') {
      if (url.searchParams.get('comp') === 'block') state.blocks.set(url.searchParams.get('blockid'), body);
      else if (url.searchParams.get('comp') === 'blocklist') {
        const ids = [...body.toString().matchAll(/<Latest>([^<]+)<\/Latest>/g)].map(match => match[1]);
        state.zip = Buffer.concat(ids.map(id => state.blocks.get(id)));
      } else throw new Error('Unexpected blob request');
      response.writeHead(201, {'Content-Length': '0', ETag: '"fixture"', 'x-ms-request-id': 'fixture'});
      response.end();return;
    }
    const method = url.pathname.split('/').pop();
    const data = JSON.parse(body);
    let result;
    if (method === 'ListArtifacts') {
      assert.equal(data.name_filter, 'ci-build-fixture');
      result = {artifacts: state.exists ? [{databaseId: '7', name: 'ci-build-fixture',
        workflowRunBackendId: 'run', workflowJobRunBackendId: 'job'}] : []};
    } else if (method === 'DeleteArtifact') {
      state.exists = false;++state.deletes;result = {ok: true, artifactId: '7'};
    } else if (method === 'CreateArtifact') {
      ++state.creates;
      if (failure === 'always-create' || (failure === 'create' && state.creates === 1)) {
        request.socket.destroy();return;
      }
      state.blocks.clear();state.zip = null;
      result = {ok: true, signedUploadUrl: `${base}/account/artifacts/file.zip?sig=fixture-secret`};
    } else if (method === 'FinalizeArtifact') {
      ++state.finalizes;state.exists = true;
      assert(state.zip?.length > 0);
      assert.equal(data.hash, 'sha256:' + crypto.createHash('sha256').update(state.zip).digest('hex'));
      assert.equal(Number(data.size), state.zip.length);
      if (failure === 'finalize' && state.finalizes === 1) {request.socket.destroy();return;}
      result = {ok: true, artifactId: '7'};
    } else throw new Error(`Unexpected RPC ${method}`);
    response.writeHead(200, {'Content-Type': 'application/json'});response.end(JSON.stringify(result));
  });
  await new Promise((resolve, reject) => {server.once('error', reject);server.listen(0, '127.0.0.1', resolve);});
  base = `http://127.0.0.1:${server.address().port}`;
  const values = {ACTIONS_RESULTS_URL: base, GITHUB_WORKSPACE: folder,
    ACTIONS_RUNTIME_TOKEN: 'e30.' + Buffer.from(JSON.stringify({scp: 'Actions.Results:run:job'})).toString('base64url') + '.fixture',
    ACTIONS_ARTIFACT_UPLOAD_CONCURRENCY: '1', ACTIONS_ARTIFACT_UPLOAD_TIMEOUT_MS: '5000'};
  const previous = Object.fromEntries(Object.keys(values).map(key => [key, process.env[key]]));
  Object.assign(process.env, values);
  t.after(async () => {
    for (const [key, value] of Object.entries(previous)) {
      if (value === undefined) delete process.env[key];else process.env[key] = value;
    }
    server.closeAllConnections();await new Promise(resolve => server.close(resolve));
    await fs.rm(folder, {recursive: true, force: true});
  });
  const outputs = {}, warnings = [];
  const core = {info: () => {}, warning: text => warnings.push(text), setOutput: (key, value) => {outputs[key] = value;}};
  const input = readInputs({UPLOAD_NAME: 'ci-build-fixture', UPLOAD_PATH: path.join(folder, 'payload.bin'),
    UPLOAD_OVERWRITE: 'true', UPLOAD_COMPRESSION: '0', UPLOAD_NO_FILES: 'error'});
  return {state, outputs, warnings, core, input, payload};
}

for (const failure of ['create', 'finalize']) {
  test(`official upload streams recover from a ${failure} connection reset`, {timeout: 30000}, async t => {
    const f = await fixture(t, failure);
    const result = await uploadArtifact(new DefaultArtifactClient(), glob, f.input, f.core, async () => {});
    assert.equal(result.id, 7);assert.equal(f.state.creates, 2);
    assert.equal(f.state.deletes, failure === 'finalize' ? 1 : 0);
    assert.equal(f.outputs['artifact-digest'], crypto.createHash('sha256').update(f.state.zip).digest('hex'));
    assert.equal(f.warnings.length, 1);
    const extracted = [];
    const parser = unzip.Parse();
    parser.on('entry', entry => {
      const chunks = [];entry.on('data', chunk => chunks.push(chunk));
      entry.on('end', () => extracted.push({name: entry.path, bytes: Buffer.concat(chunks)}));
    });
    await new Promise((resolve, reject) => {
      parser.on('error', reject);parser.on('close', resolve);
      Readable.from([f.state.zip]).pipe(parser);
    });
    assert.deepEqual(extracted, [{name: 'payload.bin', bytes: f.payload}]);
  });
}

test('official repeated CreateArtifact resets exhaust the bound without outputs', {timeout: 30000}, async t => {
  const f = await fixture(t, 'always-create');
  await assert.rejects(uploadArtifact(new DefaultArtifactClient(), glob, f.input, f.core, async () => {}), /ECONNRESET/);
  assert.equal(f.state.creates, 3);assert.equal(f.state.finalizes, 0);
  assert.deepEqual(f.outputs, {});
});

test('composite upload script runs with github-script injected parameters', {timeout: 30000}, async t => {
  const f = await fixture(t, 'create');
  const runner = await fs.mkdtemp(path.join(os.tmpdir(), 'rodin-upload-action-'));
  await fs.symlink(root, path.join(runner, 'rodin-artifact-client'), 'dir');
  const values = {RUNNER_TEMP: runner,
    UPLOAD_ACTION_PATH: path.resolve('.github/actions/upload-artifacts'),
    UPLOAD_NAME: f.input.name, UPLOAD_PATH: f.input.patterns,
    UPLOAD_OVERWRITE: 'true', UPLOAD_NO_FILES: 'error', UPLOAD_COMPRESSION: '0',
    UPLOAD_RETENTION: '0', UPLOAD_HIDDEN: 'false'};
  const previous = Object.fromEntries(Object.keys(values).map(key => [key, process.env[key]]));
  Object.assign(process.env, values);
  t.after(async () => {
    for (const [key, value] of Object.entries(previous)) {
      if (value === undefined) delete process.env[key];else process.env[key] = value;
    }
    await fs.rm(runner, {recursive: true, force: true});
  });
  const action = await fs.readFile('.github/actions/upload-artifacts/action.yml', 'utf8');
  const script = action.split('script: |\n')[1].split('\n').map(line => line.replace(/^          /, '')).join('\n');
  const AsyncFunction = Object.getPrototypeOf(async () => {}).constructor;
  const execute = new AsyncFunction('require', '__original_require__', 'github', 'octokit',
    'context', 'core', 'exec', 'glob', 'io', script);
  await execute(require, require, {}, {}, {}, f.core, {}, {}, {});
  assert.equal(f.outputs['artifact-id'], 7);
  assert.equal(f.state.creates, 2);
});
