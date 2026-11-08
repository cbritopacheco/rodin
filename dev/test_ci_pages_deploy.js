'use strict';

const assert = require('node:assert/strict');
const {test} = require('node:test');
const {deployPages} = require('../.github/actions/publish-pages/deploy.js');
const version = 'a'.repeat(40);

function fixture(statuses = ['succeed']) {
  const requests = [], output = {}, logs = [];
  let clock = 0;
  const github = {async request(route, parameters) {
    requests.push({route, parameters});
    if (route.startsWith('POST') && !route.endsWith('/cancel')) {
      return {data: {id: parameters.pages_build_version, page_url: 'https://owner.github.io/repo/'}};
    }
    if (route.endsWith('/cancel')) return {data: {}};
    const status = statuses.length > 1 ? statuses.shift() : statuses[0];
    if (status instanceof Error) throw status;
    return {data: {status}};
  }};
  const core = {async getIDToken() { return 'secret-oidc'; },
    setOutput(key, value) { output[key] = value; },
    info(message) { logs.push(message); }, warning(message) { logs.push(message); }};
  const options = {now: () => clock, sleep: async ms => { clock += ms; }, timeout: 15000};
  return {github, core, options, requests, output, logs};
}

test('distinct snapshots deploy distinct versions even with the same source commit', async () => {
  const f = fixture(['queued', 'succeed']);
  for (const [id, snapshot] of [['17', version], ['18', 'b'.repeat(40)]]) {
    assert.equal(await deployPages(f.github, f.core, 'owner/repo', id, snapshot, f.options),
      'https://owner.github.io/repo/');
  }
  const creates = f.requests.filter(x => x.route.endsWith('/deployments'));
  assert.deepEqual(creates.map(x => x.parameters.pages_build_version), [version, 'b'.repeat(40)]);
  assert.deepEqual(creates.map(x => x.parameters.artifact_id), [17, 18]);
  assert.ok(creates.every(x => x.parameters.oidc_token === 'secret-oidc'));
  assert.equal(f.output.page_url, 'https://owner.github.io/repo/');
  assert.ok(!f.logs.join('\n').includes('secret-oidc'));
});

test('permanent deployment failure is not publication success', async () => {
  const f = fixture(['deployment_perms_error']);
  await assert.rejects(deployPages(f.github, f.core, 'owner/repo', 17, version, f.options), /perms_error/);
  assert.equal(f.output.page_url, undefined);
});

test('bounded status polling cancels a stuck deployment', async () => {
  const f = fixture(['queued']);
  const listeners = process.listenerCount('SIGTERM');
  await assert.rejects(deployPages(f.github, f.core, 'owner/repo', 17, version, f.options), /Timed out/);
  assert.equal(f.requests.filter(x => x.route.endsWith('/cancel')).length, 1);
  assert.equal(f.output.page_url, undefined);
  assert.equal(process.listenerCount('SIGTERM'), listeners);
});

test('transient status error is retried, authorization error cancels and fails', async () => {
  const transient = Object.assign(new Error('Unavailable'), {status: 503});
  const f = fixture([transient, 'succeed']);
  await deployPages(f.github, f.core, 'owner/repo', 17, version, f.options);
  const g = fixture([Object.assign(new Error('Forbidden'), {status: 403})]);
  await assert.rejects(deployPages(g.github, g.core, 'owner/repo', 17, version, g.options), /Forbidden/);
  assert.equal(g.requests.filter(x => x.route.endsWith('/cancel')).length, 1);
});

test('invalid snapshot identity fails before requesting a token or deployment', async () => {
  const f = fixture();
  for (const args of [['owner/repo', 0, version], ['owner/repo', 17, 'source'],
    ['../repo', 17, version], ['owner/repo', '9007199254740993', version]]) {
    await assert.rejects(deployPages(f.github, f.core, ...args, f.options), /Expected/);
  }
  assert.equal(f.requests.length, 0);
});
