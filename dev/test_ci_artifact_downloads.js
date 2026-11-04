'use strict';

const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const os = require('node:os');
const path = require('node:path');
const {test} = require('node:test');
const {downloadArtifacts, readNames} = require('../.github/actions/download-artifacts/download.js');
const core = {info() {}, warning() {}};

async function temporary(run) {
  const folder = await fs.mkdtemp(path.join(os.tmpdir(), 'rodin-artifact-test-'));
  try { await run(folder); }
  finally { await fs.rm(folder, {recursive: true, force: true}); }
}

test('exact-name validation rejects globs, traversal, duplicates and ambiguous input', async () => {
  for (const names of [[], ['*'], ['../escape'], ['a', 'a'], ['a/b'], [1], {}]) {
    await assert.rejects(readNames({ARTIFACT_NAMES: JSON.stringify(names)}));
  }
  await assert.rejects(readNames({}));
  await assert.rejects(readNames({ARTIFACT_NAMES: '["a"]', ARTIFACT_NAMES_FILE: 'also.json'}));
  assert.deepEqual(await readNames({ARTIFACT_NAMES: '["ci-build-a"]'}), ['ci-build-a']);
  await temporary(async folder => {
    const file = path.join(folder, 'names.json');
    await fs.writeFile(file, '["ci-build-b"]');
    assert.deepEqual(await readNames({ARTIFACT_NAMES_FILE: file}), ['ci-build-b']);
  });
});

test('downloads every exact name beyond 1000 without listing or REST findBy', async () => {
  await temporary(async folder => {
    const names = Array.from({length: 1005}, (_, i) => `ci-result-${i}`);
    let gets = 0, downloads = 0;
    const client = {
      async getArtifact(name, options) {
        assert.equal(options, undefined);
        assert.equal(name, names[gets++]);
        return {artifact: {name, id: gets}};
      },
      async downloadArtifact(id, options) {
        assert.equal(options.findBy, undefined);
        assert.equal(id, ++downloads);
        await fs.writeFile(path.join(options.path, 'receipt.json'), String(id));
      },
      async listArtifacts() { assert.fail('Must never list all artifacts'); }
    };
    await downloadArtifacts(client, names, folder, core);
    assert.equal(downloads, 1005);
    assert.equal(await fs.readFile(path.join(folder, names[1004], 'receipt.json'), 'utf8'), '1005');
    assert.equal((await fs.readdir(folder)).length, 1005);
  });
});

test('failed transfer retries a resolved ID with clean staging', async () => {
  await temporary(async folder => {
    let gets = 0, attempts = 0;
    const sleeps = [];
    const client = {
      async getArtifact(name) { ++gets; return {artifact: {name, id: 7}}; },
      async downloadArtifact(id, {path: staging}) {
        assert.equal(id, 7);
        assert.deepEqual(await fs.readdir(staging), []);
        if (++attempts === 1) {
          await fs.writeFile(path.join(staging, 'partial'), 'bad');
          throw new Error('Interrupted transfer');
        }
        await fs.writeFile(path.join(staging, 'receipt.json'), 'complete');
      }
    };
    await downloadArtifacts(client, ['ci-target-a'], folder, core, async delay => sleeps.push(delay));
    assert.equal(gets, 1);
    assert.equal(attempts, 2);
    assert.deepEqual(sleeps, [2000]);
    assert.deepEqual(await fs.readdir(path.join(folder, 'ci-target-a')), ['receipt.json']);
  });
});

test('missing artifact fails without silently accepting an empty selection', async () => {
  await temporary(async folder => {
    const client = {async getArtifact() { throw new Error('Artifact not found'); }};
    await assert.rejects(downloadArtifacts(client, ['ci-target-missing'], folder, core), /not found/);
    assert.deepEqual(await fs.readdir(folder), []);
  });
});

test('empty or permanently failed downloads exhaust retries and remain failures', async () => {
  for (const fail of [false, true]) await temporary(async folder => {
    let attempts = 0;
    const client = {
      async getArtifact(name) { return {artifact: {name, id: 8}}; },
      async downloadArtifact() { ++attempts; if (fail) throw new Error('Unavailable'); }
    };
    await assert.rejects(downloadArtifacts(client, ['ci-target-a'], folder, core, async () => {}));
    assert.equal(attempts, 3);
    assert.deepEqual(await fs.readdir(folder), []);
  });
});

test('incorrect resolved name fails before download', async () => {
  await temporary(async folder => {
    const client = {async getArtifact() { return {artifact: {name: 'other', id: 1}}; }};
    await assert.rejects(downloadArtifacts(client, ['expected'], folder, core), /wrong name/);
  });
});

test('digest mismatch is retried and never published as valid input', async () => {
  await temporary(async folder => {
    let attempts = 0;
    const client = {
      async getArtifact(name) { return {artifact: {name, id: 9, digest: 'sha256:expected'}}; },
      async downloadArtifact(id, {path: staging, expectedHash}) {
        assert.equal(expectedHash, 'sha256:expected');
        assert.deepEqual(await fs.readdir(staging), []);
        await fs.writeFile(path.join(staging, 'payload'), 'corrupted');
        ++attempts;
        return {digestMismatch: true};
      }
    };
    await assert.rejects(downloadArtifacts(client, ['benchmark-group-a'], folder, core,
      async () => {}), /digest mismatch/);
    assert.equal(attempts, 3);
    assert.deepEqual(await fs.readdir(folder), []);
  });
});

test('lookup resets recover within three attempts and preserve the resolved ID', async () => {
  await temporary(async folder => {
    let gets = 0, downloads = 0;
    const delays = [];
    const client = {
      async getArtifact(name) {
        if (++gets < 3) throw new Error('Failed to ListArtifacts: Unable to make request: ECONNRESET');
        return {artifact: {name, id: 42}};
      },
      async downloadArtifact(id, {path: staging}) {
        ++downloads;
        assert.equal(id, 42);
        await fs.writeFile(path.join(staging, 'index.json'), 'verified');
      }
    };
    await downloadArtifacts(client, ['ci-build-index-a'], folder, core, async ms => delays.push(ms));
    assert.equal(gets, 3);
    assert.equal(downloads, 1);
    assert.deepEqual(delays, [2000, 4000]);
    assert.equal(await fs.readFile(path.join(folder, 'ci-build-index-a/index.json'), 'utf8'), 'verified');
  });
});

test('lookup exhaustion and permanent failures never publish input', async () => {
  for (const [error, expected] of [
    [new Error('ECONNRESET'), 3],
    [Object.assign(new Error('service unavailable'), {statusCode: 503}), 3],
    [Object.assign(new Error('unauthorized'), {statusCode: 403}), 1],
    [Object.assign(new Error('missing'), {name: 'ArtifactNotFoundError'}), 1],
    [new Error('unknown failure'), 1],
    [Object.assign(new Error('ECONNRESET'), {name: 'AbortError'}), 1]
  ]) await temporary(async folder => {
    let attempts = 0;
    const client = {
      async getArtifact() { ++attempts; throw error; },
      async downloadArtifact() { assert.fail('Unresolved artifact cannot be downloaded'); }
    };
    await assert.rejects(downloadArtifacts(client, ['ci-build-index-a'], folder, core, async () => {}));
    assert.equal(attempts, expected);
    assert.deepEqual(await fs.readdir(folder), []);
  });
});

test('single-artifact merge preserves workspace and publishes only verified bytes', async () => {
  await temporary(async folder => {
    await fs.writeFile(path.join(folder, 'source.cpp'), 'keep');
    await fs.writeFile(path.join(folder, 'ci-build-index.json'), 'old');
    let gets = 0, transfers = 0;
    const client = {
      async getArtifact(name) {
        if (++gets === 1) throw Object.assign(new Error('reset'), {code: 'ECONNRESET'});
        return {artifact: {name, id: 51, digest: 'sha256:verified'}};
      },
      async downloadArtifact(id, {path: staging, expectedHash}) {
        assert.equal(expectedHash, 'sha256:verified');
        assert.equal(await fs.readFile(path.join(folder, 'ci-build-index.json'), 'utf8'), 'old');
        assert.deepEqual(await fs.readdir(staging), []);
        if (++transfers === 1) {
          await fs.writeFile(path.join(staging, 'partial'), 'bad');
          throw new Error('Interrupted transfer');
        }
        await fs.writeFile(path.join(staging, 'ci-build-index.json'), 'new');
        await fs.mkdir(path.join(staging, 'nested'));
        await fs.writeFile(path.join(staging, 'nested/receipt'), 'valid');
      }
    };
    await downloadArtifacts(client, ['ci-build-index-a'], folder, core, async () => {}, true);
    assert.equal(gets, 2);
    assert.equal(transfers, 2);
    assert.equal(await fs.readFile(path.join(folder, 'source.cpp'), 'utf8'), 'keep');
    assert.equal(await fs.readFile(path.join(folder, 'ci-build-index.json'), 'utf8'), 'new');
    assert.equal(await fs.readFile(path.join(folder, 'nested/receipt'), 'utf8'), 'valid');
    assert.deepEqual((await fs.readdir(folder)).sort(), ['ci-build-index.json', 'nested', 'source.cpp']);
    await assert.rejects(downloadArtifacts(client, ['one', 'two'], folder, core, async () => {}, true), /exactly one/);
  });
});

test('failed merged download leaves existing workspace files unchanged', async () => {
  await temporary(async folder => {
    await fs.writeFile(path.join(folder, 'index.json'), 'old');
    const client = {
      async getArtifact(name) { return {artifact: {name, id: 52}}; },
      async downloadArtifact(id, {path: staging}) {
        await fs.writeFile(path.join(staging, 'index.json'), 'bad');
        return {digestMismatch: true};
      }
    };
    await assert.rejects(downloadArtifacts(client, ['ci-build-index-a'], folder, core, async () => {}, true), /digest mismatch/);
    assert.equal(await fs.readFile(path.join(folder, 'index.json'), 'utf8'), 'old');
    assert.deepEqual(await fs.readdir(folder), ['index.json']);
  });
});
