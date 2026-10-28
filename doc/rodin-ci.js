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
  const workflows = ['Tests', 'Benchmarks', 'Build', 'Style', 'CodeQL', 'Labeler', 'Labeler Rules', 'Copilot setup steps', 'Documentation', 'Publish documentation manually'];
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
      status.textContent = selected + ' · ' + sha + ' · Checked ' + new Date().toLocaleString();
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
        if (run?.head_sha === sha && run.conclusion === 'success') result.className = 'm-text m-success';
        if (run?.conclusion === 'failure') result.className = 'm-text m-danger';
        row.insertCell().textContent = run ? run.head_sha.slice(0, 10) : '—';
        const title = run && 'CI report · ' + name + ' · ' + selected + ' · ' + run.head_sha + ' · attempt ' + run.run_attempt;
        const report = title && reports.workflow_runs.find(item => item.display_title === title && item.conclusion === 'success');
        row.insertCell().append(report ? link('Open report', report.html_url) : link('Run summary', run?.html_url));
      }
      output.append(table);
      const note = document.createElement('p');
      note.textContent = 'Open report links to a completed generated report when available; Run summary opens GitHub’s native run page. Missing recent runs are unknown, not successful. Full coverage and campaign completeness are confirmed in the reports.';
      output.append(note);
    } catch (error) {
      status.textContent = 'CI status unavailable: ' + error.message + '. Use the workflow links below.';
      for (const name of workflows) {
        const p = document.createElement('p'); p.append(link(name, actions + '?query=' + encodeURIComponent('workflow:"' + name + '"'))); output.append(p);
      }
    } finally {refresh.disabled = branch.disabled = false;}
  }
  refresh.addEventListener('click', load); branch.addEventListener('change', load); load();
})();
