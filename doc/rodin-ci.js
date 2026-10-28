/* Public CI index. Workflow results and report links are never inferred across commits. */
(() => {
  if (!location.pathname.endsWith('/ci-index.html')) return;
  const section = document.getElementById('ci-status');
  if (!section) return;
  const root = document.createElement('div');
  section.append(root);
  const branch = document.createElement('select');
  branch.setAttribute('aria-label', 'CI branch');
  for (const name of ['develop', 'master']) branch.add(new Option(name, name));
  const refresh = document.createElement('button');
  refresh.type = 'button'; refresh.textContent = 'Refresh'; refresh.className = 'm-button m-default';
  const status = document.createElement('p'); status.setAttribute('role', 'status');
  const output = document.createElement('div');
  root.append(branch, ' ', refresh, status, output);
  const api = 'https://api.github.com/repos/cbritopacheco/rodin';
  const actions = 'https://github.com/cbritopacheco/rodin/actions';
  const workflows = ['Tests', 'Benchmarks', 'Build', 'Style', 'CodeQL', 'Labeler', 'Labeler Rules', 'Copilot setup steps', 'Documentation'];
  const link = (text, url) => {
    const node = document.createElement('a'); node.textContent = text;
    // API metadata is untrusted: link only to this repository's Actions pages.
    node.href = typeof url === 'string' && url.startsWith(actions + '/') ? url : actions;
    return node;
  };
  async function read(path) {
    const response = await fetch(api + path, {headers: {Accept: 'application/vnd.github+json'}});
    if (!response.ok) throw new Error('GitHub API returned HTTP ' + response.status);
    return response.json();
  }
  async function load() {
    refresh.disabled = branch.disabled = true;
    const selected = branch.value;
    status.textContent = 'Loading ' + selected + ' CI…'; output.replaceChildren();
    try {
      const [head, runs, reports] = await Promise.all([
        read('/branches/' + selected),
        read('/actions/runs?branch=' + selected + '&per_page=100'),
        read('/actions/workflows/RetryInterruptedJobs.yml/runs?status=completed&per_page=100').catch(() => ({workflow_runs: []}))
      ]);
      const sha = head.commit.sha;
      status.textContent = selected + ' · ';
      const headLink = document.createElement('a');
      headLink.textContent = sha; headLink.href = 'https://github.com/cbritopacheco/rodin/commit/' + sha;
      status.append(headLink, ' · Checked ' + new Date().toLocaleString());
      const table = document.createElement('table'); table.className = 'm-table';
      const header = table.createTHead().insertRow();
      for (const text of ['Workflow', 'Status', 'Commit', 'Report']) {
        const th = document.createElement('th'); th.textContent = text; header.append(th);
      }
      const body = table.createTBody();
      for (const name of workflows) {
        const run = runs.workflow_runs.find(item => item.name === name);
        const row = body.insertRow();
        row.insertCell().append(link(name, run?.html_url));
        const result = row.insertCell();
        result.textContent = run ? (run.conclusion || run.status) : 'No recent run';
        if (run && run.head_sha !== sha) result.append(' · older commit');
        // Green is reserved for a successful result on the current branch head.
        const outcome = run?.conclusion || run?.status;
        result.className = 'rodin-ci-status ' + (
          run?.head_sha !== sha ? 'rodin-ci-muted' :
          outcome === 'success' ? 'rodin-ci-success' :
          ['failure', 'timed_out', 'startup_failure', 'action_required'].includes(outcome) ? 'rodin-ci-danger' :
          ['queued', 'in_progress', 'waiting', 'pending', 'requested'].includes(outcome) ? 'rodin-ci-warning' : 'rodin-ci-muted');
        const commit = row.insertCell();
        if (run && /^[a-f0-9]{40}$/i.test(run.head_sha)) {
          const commitLink = document.createElement('a');
          commitLink.textContent = run.head_sha.slice(0, 10);
          commitLink.href = 'https://github.com/cbritopacheco/rodin/commit/' + run.head_sha;
          commit.append(commitLink);
        } else commit.textContent = '—';
        const title = run && 'CI report · ' + name + ' · ' + selected + ' · ' + run.head_sha + ' · attempt ' + run.run_attempt;
        const report = title && reports.workflow_runs.find(item => item.display_title === title && item.conclusion === 'success');
        row.insertCell().append(report ? link('Open report', report.html_url) : link('Summary', run?.html_url));
      }
      const actionsLink = document.createElement('p');
      actionsLink.style.textAlign = 'center';
      actionsLink.append(link('Open GitHub Actions', actions));
      output.append(table, actionsLink);
    } catch (error) {
      status.textContent = 'CI status is currently unavailable. Open a workflow below to check its status.';
      for (const name of workflows) {
        const p = document.createElement('p'); p.append(link(name, actions + '?query=' + encodeURIComponent('workflow:"' + name + '"'))); output.append(p);
      }
    } finally {refresh.disabled = branch.disabled = false;}
  }
  refresh.addEventListener('click', load); branch.addEventListener('change', load); load();
})();
