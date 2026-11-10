'use strict';

// Use the published gh-pages commit to distinguish master/develop/benchmark
// snapshots produced by a shared source commit.
async function deployPages(github, core, repository, artifactId, version,
  {sleep = ms => new Promise(resolve => setTimeout(resolve, ms)),
    now = Date.now, timeout = 600000} = {}) {
  if (!/^[\w-]+\/[\w.-]+$/.test(repository) ||
      !/^[a-f0-9]{40}$/.test(version) || !/^[1-9][0-9]*$/.test(String(artifactId)) ||
      !Number.isSafeInteger(Number(artifactId))) {
    throw new Error('Expected a repository, artifact ID and exact gh-pages commit');
  }
  const [owner, repo] = repository.split('/');
  const {data: deployment} = await github.request('POST /repos/{owner}/{repo}/pages/deployments', {
    owner, repo, artifact_id: Number(artifactId), pages_build_version: version,
    oidc_token: await core.getIDToken()
  });
  const deploymentId = deployment.id || version;
  let pending = true;
  async function cancel() {
    if (!pending) return;
    await github.request('POST /repos/{owner}/{repo}/pages/deployments/{deploymentId}/cancel',
      {owner, repo, deploymentId});
    pending = false;
  }
  const onSignal = () => cancel().catch(error => core.warning(error.message))
    .finally(() => process.exit(1));
  process.once('SIGINT', onSignal);
  process.once('SIGTERM', onSignal);
  try {
    const deadline = now() + timeout;
    let errors = 0;
    while (now() < deadline) {
      let status;
      try {
        const response = await github.request('GET /repos/{owner}/{repo}/pages/deployments/{deploymentId}',
          {owner, repo, deploymentId});
        status = response.data.status;
        errors = 0;
      } catch (error) {
        if (![404, 429, 500, 502, 503, 504].includes(error.status) || ++errors >= 10) throw error;
        core.warning(`Temporary Pages status failure (${error.status})`);
      }
      if (status === 'succeed') {
        pending = false;
        core.setOutput('page_url', deployment.page_url);
        core.info(`Pages deployed snapshot ${version}`);
        return deployment.page_url;
      }
      if (['deployment_failed', 'deployment_perms_error', 'deployment_content_failed',
        'deployment_cancelled', 'deployment_lost'].includes(status)) {
        pending = false;
        throw new Error(`Pages deployment failed: ${status}`);
      }
      await sleep(5000);
    }
    throw new Error(`Timed out deploying Pages snapshot ${version}`);
  } catch (error) {
    try { await cancel(); }
    catch (cleanup) { core.warning(`Could not cancel Pages deployment: ${cleanup.message}`); }
    throw error;
  } finally {
    process.removeListener('SIGINT', onSignal);
    process.removeListener('SIGTERM', onSignal);
  }
}

module.exports = {deployPages};
