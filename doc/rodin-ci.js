/* Public CI index. Workflow results and report links are never inferred across commits. */
(() => {
  if (!location.pathname.endsWith('/ci-index.html')) return;
  const section = document.getElementById('ci-status');
  if (!section) return;
  const api = 'https://api.github.com/repos/cbritopacheco/rodin';
  const actions = 'https://github.com/cbritopacheco/rodin/actions';
  const workflows = ['Tests', 'Benchmarks', 'Build', 'Style', 'CodeQL', 'Labeler', 'Preflight', 'Copilot setup steps', 'Documentation'];
  const badgeStates = {passing: 'success', failing: 'failure', running: 'in_progress',
    queued: 'queued', cancelled: 'cancelled', skipped: 'skipped', unknown: 'unknown'};
  const statusNames = {success: 'passing', failure: 'failing', in_progress: 'running'};
  const link = (text, url) => {
    const node = document.createElement('a'); node.textContent = text;
    // API metadata is untrusted: link only to this repository's Actions pages.
    try {
      const target = new URL(url);
      node.href = target.origin === 'https://github.com' && !target.username && !target.password &&
        (target.pathname === '/cbritopacheco/rodin/actions' || target.pathname.startsWith('/cbritopacheco/rodin/actions/'))
        ? target.href : actions;
    } catch {node.href = actions;}
    return node;
  };
  async function read(path) {
    const response = await fetch(api + path, {headers: {Accept: 'application/vnd.github+json'}});
    if (!response.ok) throw new Error('GitHub API returned HTTP ' + response.status);
    return response.json();
  }
  async function badge(branch, name) {
    try {
      const response = await fetch('https://raw.githubusercontent.com/cbritopacheco/rodin/ci-badges/' +
        branch + '/' + name.toLowerCase() + '.json', {cache: 'no-store'});
      if (!response.ok) return null;
      const value = await response.json();
      return value.schemaVersion === 1 && Object.hasOwn(badgeStates, value.message) ? value : null;
    } catch {return null;}
  }
  for (const selected of ['master', 'develop']) {
    const root = document.createElement('section');
    root.className = 'rodin-ci-branch';
    root.id = 'ci-' + selected;
    const heading = document.createElement('h3'); heading.textContent = selected;
    const refresh = document.createElement('button');
    refresh.type = 'button'; refresh.textContent = 'Refresh'; refresh.className = 'm-button m-default';
    refresh.setAttribute('aria-label', 'Refresh ' + selected + ' CI');
    const status = document.createElement('p'); status.setAttribute('role', 'status');
    status.className = 'rodin-ci-checked';
    const output = document.createElement('div');
    root.append(heading, refresh, status, output);
    section.append(root);
    async function load() {
      refresh.disabled = true;
      status.textContent = 'Loading ' + selected + ' CI…'; output.replaceChildren();
      try {
        const [head, runs, reports, badges] = await Promise.all([
          read('/branches/' + selected),
          read('/actions/runs?branch=' + selected + '&per_page=100'),
          read('/actions/workflows/RetryInterruptedJobs.yml/runs?status=completed&per_page=100').catch(() => ({workflow_runs: []})),
          Promise.all(['Tests', 'Benchmarks', 'Build', 'Documentation'].map(async name => [name, await badge(selected, name)]))
        ]);
        const badgeByName = Object.fromEntries(badges);
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
          row.insertCell().append(link(name === 'Labeler' ? 'Labeler (pull requests)' : name, run?.html_url));
          const result = row.insertCell();
          const published = badgeByName[name];
          // The publisher resolves job activity and verifies documentation deployment.
          // An older badge never establishes status for the displayed branch head.
          const sameRun = name !== 'Documentation' || (published?.run === String(run?.id) &&
            published?.attempt === run?.run_attempt);
          const outcome = run?.head_sha === sha && published?.commit === sha && sameRun
            ? badgeStates[published.message] : run?.conclusion || run?.status;
          result.textContent = run ? (statusNames[outcome] || outcome.replaceAll('_', ' ')) : 'No recent run';
          if (run && run.head_sha !== sha) result.append(' · older commit');
          // Green is reserved for a successful result on the current branch head.
          result.className = 'rodin-ci-status ' + (
            run?.head_sha !== sha ? 'rodin-ci-muted' :
            outcome === 'success' ? 'rodin-ci-success' :
            ['failure', 'timed_out', 'startup_failure', 'action_required'].includes(outcome) ? 'rodin-ci-danger' :
            outcome === 'in_progress' ? 'rodin-ci-warning' :
            ['queued', 'waiting', 'pending', 'requested'].includes(outcome) ? 'rodin-ci-queued' :
            ['cancelled', 'skipped'].includes(outcome) ? 'rodin-ci-cancelled' : 'rodin-ci-muted');
          const commit = row.insertCell();
          if (run && /^[a-f0-9]{40}$/i.test(run.head_sha)) {
            const commitLink = document.createElement('a');
            commitLink.textContent = run.head_sha.slice(0, 10);
            commitLink.href = 'https://github.com/cbritopacheco/rodin/commit/' + run.head_sha;
            commit.append(commitLink);
          } else commit.textContent = '—';
          const title = run && 'CI report · ' + name + ' · ' + selected + ' · ' + run.head_sha + ' · attempt ' + run.run_attempt;
          const report = title && reports.workflow_runs.find(item => item.display_title === title && item.conclusion === 'success');
          const links = row.insertCell();
          if (run) links.append(link('Workflow run', run.html_url));
          if (report) {
            // Reports have one reporting job; its check-run ID identifies the
            // generated summary on GitHub's workflow overview.
            const jobs = await read('/actions/runs/' + report.id + '/jobs?filter=latest&per_page=100')
              .catch(() => ({jobs: []}));
            const job = jobs.jobs.find(item => item.name === 'Inspect interruptions and request one bounded retry');
            if (job && Number.isSafeInteger(job.id)) {
              links.append(' · ', link('Report', report.html_url + '/attempts/' + report.run_attempt + '#summary-' + job.id));
            }
          }
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
      } finally {refresh.disabled = false;}
    }
    refresh.addEventListener('click', load); load();
  }
})();
