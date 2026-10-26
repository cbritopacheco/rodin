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
  const server = http.createServer((request, response) => {
    response.writeHead(200, {'Content-Type': 'application/zip', 'Content-Length': zip.length});
    // Close a valid ZIP stream before its advertised payload has arrived.
    response.end(interrupted ? zip.subarray(0, 1024) : zip);
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
