'use strict';

// Exercise the pinned official client's real streams, not a mocked download.
const assert = require('node:assert/strict');
const {execFile} = require('node:child_process');
const crypto = require('node:crypto');
const fs = require('node:fs/promises');
const http = require('node:http');
const os = require('node:os');
const path = require('node:path');
const {test} = require('node:test');

const clientRoot = process.env.ARTIFACT_CLIENT_ROOT ||
  path.resolve('.github/actions/download-artifacts');
const clientModule = path.join(clientRoot,
  'node_modules/@actions/artifact/lib/internal/download/download-artifact.js');
const {DefaultArtifactClient} = require(path.join(clientRoot, 'node_modules/@actions/artifact'));
const {streamExtractExternal} = require(clientModule);
const archiver = require(path.join(clientRoot, 'node_modules/archiver'));

async function fixture(t, interrupted) {
  const folder = await fs.mkdtemp(path.join(os.tmpdir(), 'rodin-artifact-stream-'));
  const payload = Buffer.alloc(65536, 17);
  const archive = archiver('zip', {store: true});
  const chunks = [];
  const finished = new Promise((resolve, reject) => {
    archive.on('data', chunk => chunks.push(chunk));
    archive.on('end', () => resolve(Buffer.concat(chunks)));
    archive.on('error', reject);
  });
  archive.append(payload, {name: 'payload.bin'});
  await archive.finalize();
  const zip = await finished;
  let requests = 0;
  const server = http.createServer((request, response) => {
    response.writeHead(200, {'Content-Type': 'application/zip', 'Content-Length': zip.length});
    // Close a valid ZIP stream before its advertised payload has arrived.
    const partial = typeof interrupted === 'function' ? interrupted(++requests) : interrupted;
    response.end(partial ? zip.subarray(0, 1024) : zip);
  });
  t.after(async () => {
    server.closeAllConnections();
    await new Promise(resolve => server.close(resolve));
    await fs.rm(folder, {recursive: true, force: true});
  });
  await new Promise((resolve, reject) => {
    server.once('error', reject);
    server.listen(0, '127.0.0.1', resolve);
  });
  return {folder, payload, zip, url: `http://127.0.0.1:${server.address().port}/artifact.zip`};
}

test('official client completes and hashes a real artifact transfer', async t => {
  assert.equal(typeof new DefaultArtifactClient().getArtifact, 'function');
  const {folder, payload, zip, url} = await fixture(t, false);
  const result = await streamExtractExternal(url, folder, {timeout: 5000});
  assert.equal(result.sha256Digest, 'sha256:' + crypto.createHash('sha256').update(zip).digest('hex'));
  assert.deepEqual(await fs.readFile(path.join(folder, 'payload.bin')), payload);
});

test('interrupted ZIP rejects instead of exiting successfully with an unresolved promise', async t => {
  const {folder, url} = await fixture(t, true);
  const script = `
    const original = global.setTimeout;
    // Keep the old client's hardcoded 30-second watchdog short in this fixture.
    global.setTimeout = (fn, ms, ...args) => original(fn, ms === 30000 ? 50 : ms, ...args);
    const {streamExtractExternal} = require(process.argv[1]);
    let settled = false;
    process.on('exit', () => console.log('settled=' + settled));
    streamExtractExternal(process.argv[2], process.argv[3], {timeout: 50}).then(() => {
      settled = true;
    }).catch(error => {
      settled = true;
      console.error(error.message);
      process.exitCode = 1;
    });`;
  const result = await new Promise(resolve => {
    execFile(process.execPath, ['-e', script, clientModule, url, folder],
      {timeout: 10000}, (error, stdout, stderr) => resolve({error, stdout, stderr}));
  });
  assert.equal(result.error && result.error.code, 1, result.stderr);
  assert.match(result.stdout, /settled=true/);
  assert.match(result.stderr, /Blob storage chunk did not respond/);
});


test('mandatory download recovers a real interrupted ZIP with clean staging', async t => {
  const {downloadArtifacts} = require('../.github/actions/download-artifacts/download.js');
  const {folder, payload, zip, url} = await fixture(t, request => request === 1);
  let attempts = 0;
  const digest = 'sha256:' + crypto.createHash('sha256').update(zip).digest('hex');
  const client = {
    async getArtifact(name) { return {artifact: {name, id: 17, digest}}; },
    async downloadArtifact(id, options) {
      assert.equal(id, 17);
      assert.equal(options.expectedHash, digest);
      assert.deepEqual(await fs.readdir(options.path), []);
      ++attempts;
      const result = await streamExtractExternal(url, options.path, {timeout: 50});
      return {digestMismatch: result.sha256Digest !== digest};
    }
  };
  await downloadArtifacts(client, ['benchmark-group-a'], folder,
    {info() {}, warning() {}}, async () => {});
  assert.equal(attempts, 2);
  assert.deepEqual(await fs.readFile(path.join(folder, 'benchmark-group-a/payload.bin')), payload);
  assert.deepEqual(await fs.readdir(folder), ['benchmark-group-a']);
});

test('composite download retries real lookup reset and merges the verified artifact', {timeout: 30000}, async t => {
  const {folder, payload, zip, url} = await fixture(t, false);
  const runner = await fs.mkdtemp(path.join(os.tmpdir(), 'rodin-download-action-'));
  await fs.symlink(clientRoot, path.join(runner, 'rodin-artifact-client'), 'dir');
  await fs.writeFile(path.join(folder, 'source.cpp'), 'keep');
  const digest = 'sha256:' + crypto.createHash('sha256').update(zip).digest('hex');
  let lookups = 0, signed = 0;
  const server = http.createServer(async (request, response) => {
    const chunks = [];
    for await (const chunk of request) chunks.push(chunk);
    const data = JSON.parse(Buffer.concat(chunks));
    const method = request.url.split('/').pop();
    let result;
    if (method === 'ListArtifacts') {
      ++lookups;
      if (lookups === 1) {request.socket.destroy();return;}
      assert(data.name_filter === 'ci-build-index-fixture' || data.id_filter === '17');
      result = {artifacts: [{databaseId: '17', name: 'ci-build-index-fixture', size: String(zip.length),
        workflowRunBackendId: 'run', workflowJobRunBackendId: 'job', digest}]};
    } else if (method === 'GetSignedArtifactURL') {
      ++signed;
      assert.equal(data.name, 'ci-build-index-fixture');
      result = {signedUrl: url};
    } else throw new Error('Unexpected RPC: ' + method);
    response.writeHead(200, {'Content-Type': 'application/json'});
    response.end(JSON.stringify(result));
  });
  await new Promise((resolve, reject) => {server.once('error', reject);server.listen(0, '127.0.0.1', resolve);});
  const values = {RUNNER_TEMP: runner, GITHUB_WORKSPACE: folder,
    ACTIONS_RESULTS_URL: `http://127.0.0.1:${server.address().port}`,
    ACTIONS_RUNTIME_TOKEN: 'e30.' + Buffer.from(JSON.stringify({scp: 'Actions.Results:run:job'})).toString('base64url') + '.fixture',
    ARTIFACT_ACTION_PATH: path.resolve('.github/actions/download-artifacts'),
    ARTIFACT_NAMES: '["ci-build-index-fixture"]', ARTIFACT_NAMES_FILE: '',
    ARTIFACT_DESTINATION: folder, ARTIFACT_MERGE_SINGLE: 'true'};
  const previous = Object.fromEntries(Object.keys(values).map(key => [key, process.env[key]]));
  Object.assign(process.env, values);
  t.after(async () => {
    for (const [key, value] of Object.entries(previous)) {
      if (value === undefined) delete process.env[key];else process.env[key] = value;
    }
    server.closeAllConnections();await new Promise(resolve => server.close(resolve));
    await fs.rm(runner, {recursive: true, force: true});
  });
  const action = await fs.readFile('.github/actions/download-artifacts/action.yml', 'utf8');
  const script = action.split('script: |\n')[1].split('\n').map(line => line.replace(/^          /, '')).join('\n');
  const AsyncFunction = Object.getPrototypeOf(async () => {}).constructor;
  const execute = new AsyncFunction('require', '__original_require__', 'github', 'octokit',
    'context', 'core', 'exec', 'glob', 'io', script);
  const warnings = [];
  await execute(require, require, {}, {}, {}, {info() {}, warning: message => warnings.push(message)}, {}, {}, {});
  assert.equal(lookups, 3);
  assert.equal(signed, 1);
  assert.equal(warnings.length, 1);
  assert.deepEqual(await fs.readFile(path.join(folder, 'payload.bin')), payload);
  assert.equal(await fs.readFile(path.join(folder, 'source.cpp'), 'utf8'), 'keep');
  assert.deepEqual((await fs.readdir(folder)).sort(), ['payload.bin', 'source.cpp']);
});
