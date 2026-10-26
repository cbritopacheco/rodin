/* Per-component histories are published only after complete-campaign validation. */
(() => {
  const root = document.querySelector('[data-benchmark-history][data-history-index]');
  if (!root) return;
  let chart, history, selectedCase, samples, generation = 0;
  let allCases = [], configurationURL, filterTimer;
  const requested = new URLSearchParams(location.search);
  let requestedCase = requested.get('case');
  const factors = { ns: 1, us: 1000, ms: 1000000, s: 1000000000 };
  const create = (tag, text) => { const node = document.createElement(tag); if (text) node.textContent = text; return node; };
  const read = async url => {
    const response = await fetch(url);
    if (!response.ok) throw new Error('Timing history could not be loaded.');
    const data = await response.json();
    if (data.schema !== 1) throw new Error('Unsupported timing history format.');
    return data;
  };
  const render = () => {
    if (!history || !selectedCase || !window.Chart) return;
    const [name, metric] = selectedCase.split('\u0000');
    const entries = history.entries;
    samples = entries.map(entry => entry.results.find(result => result.name === name && result.metric === metric));
    if (chart) chart.destroy();
    const colors = getComputedStyle(document.documentElement);
    chart = new Chart(canvas, {
      type: 'line',
      data: { labels: entries.map(entry => entry.commit.slice(0, 10)), datasets: [{
        label: name, data: samples.map(value => value ? value.value * factors[value.unit] : null),
        borderColor: colors.getPropertyValue('--link-color').trim(), backgroundColor: 'transparent',
        pointRadius: 3, lineTension: 0, spanGaps: false
      }] },
      options: { responsive: true, maintainAspectRatio: false, legend: { display: false },
        scales: { xAxes: [{ ticks: { fontColor: colors.getPropertyValue('--color').trim(), maxTicksLimit: 8 },
                            scaleLabel: { display: true, labelString: 'Commit', fontColor: colors.getPropertyValue('--color').trim() } }],
                  yAxes: [{ ticks: { beginAtZero: true, fontColor: colors.getPropertyValue('--color').trim() },
                            scaleLabel: { display: true, labelString: metric + ' · ns / ' + (name.includes('/build/') ? 'construction' : history.format === 'csv' ? 'query' : 'iteration'), fontColor: colors.getPropertyValue('--color').trim() } }] },
        tooltips: { callbacks: { afterLabel: item => {
          const entry = entries[item.index], sample = samples[item.index];
          return [entry.date, 'Samples: ' + sample.samples.length, 'Run: ' + entry.run];
        } } },
        onClick: (_event, active) => { if (active.length) window.open('https://github.com/cbritopacheco/rodin/commit/' + entries[active[0]._index].commit, '_blank', 'noopener'); }
      }
    });
  };
  const configurationLabel = create('label', 'Configuration');
  const configurations = create('select');configurations.setAttribute('aria-label', 'Benchmark configuration');configurationLabel.append(configurations);
  const caseLabel = create('label', 'Benchmark case');
  const cases = create('select');cases.setAttribute('aria-label', 'Benchmark case');caseLabel.append(cases);
  const filterLabel = create('label', 'Filter cases');
  const filter = create('input');filter.type = 'search';filter.setAttribute('aria-label', 'Filter benchmark cases');
  filter.disabled = true;filter.placeholder = 'Operation, geometry, degree…';filter.value = requested.get('filter') || '';filterLabel.append(filter);
  const controls = create('div');controls.className = 'benchmark-chart-controls';controls.append(configurationLabel, filterLabel, caseLabel);
  const chartBox = create('div');chartBox.className = 'benchmark-chart';
  const canvas = create('canvas');canvas.setAttribute('role', 'img');canvas.setAttribute('aria-label', 'Timing history for the selected benchmark case');chartBox.append(canvas);
  const status = create('p');status.className = 'benchmark-chart-status';status.setAttribute('aria-live', 'polite');
  const fail = error => { status.textContent = error.message;status.setAttribute('role', 'alert'); };
  const loadCase = async () => {
    const request = ++generation;
    const option = cases.selectedOptions[0];
    const loaded = await read(option.dataset.url);
    if (request !== generation) return;
    history = loaded;
    const url = new URL(location.href);
    url.searchParams.set('configuration', configurations.value);
    url.searchParams.set('case', option.value);
    if (filter.value.trim()) url.searchParams.set('filter', filter.value.trim());
    else url.searchParams.delete('filter');
    window.history.replaceState(window.history.state, '', url);
    selectedCase = history.name + '\u0000' + history.metric;
    status.textContent = `Showing up to ${history.window} recent commits. Repetitions, counters and context remain in the published series.`;
    render();
  };
  const filterCases = async (preferred = cases.value) => {
    ++generation;
    const terms = filter.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
    const matching = allCases.filter(item => terms.every(term => (item.name + ' ' + item.metric).toLowerCase().includes(term)));
    cases.replaceChildren();cases.disabled = !matching.length;
    for (const item of matching) {
      const option = create('option', item.name + ' · ' + item.metric);
      option.value = item.id;option.dataset.url = new URL(item.url, configurationURL).href;cases.append(option);
    }
    if (!matching.length) {
      history = undefined;
      if (chart) chart.destroy();chart = undefined;
      status.textContent = 'No matching benchmark cases.';
      return;
    }
    if (matching.some(item => item.id === preferred)) cases.value = preferred;
    await loadCase();
  };
  const load = async () => {
    const request = ++generation;
    filter.disabled = true;cases.disabled = true;
    const option = configurations.selectedOptions[0];
    const url = new URL(option.dataset.url);
    const index = await read(url);
    if (request !== generation) return;
    if (!index.cases.length) throw new Error('No benchmark case is available.');
    configurationURL = url;
    filter.disabled = false;
    allCases = [...index.cases].sort((a,b) => a.name.localeCompare(b.name) || a.metric.localeCompare(b.metric));
    const preferred = requestedCase || cases.value;
    requestedCase = null;
    await filterCases(preferred);
  };
  filter.addEventListener('input', () => {
    ++generation;
    clearTimeout(filterTimer);
    filterTimer = setTimeout(() => filterCases().catch(fail), 150);
  });
  cases.addEventListener('change', () => loadCase().catch(fail));
  configurations.addEventListener('change', () => load().catch(fail));
  new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  matchMedia('(prefers-color-scheme: light)').addEventListener('change', render);
  window.addEventListener('load', () => { if (!window.Chart) fail(new Error('The chart renderer could not be loaded.')); else render(); }, { once: true });
  (async () => {
    const url = new URL(root.dataset.historyIndex, location.href);
    const index = await read(url);
    if (!index.configurations.length) throw new Error('No timing configuration is available.');
    root.replaceChildren(controls, chartBox, status);
    for (const configuration of index.configurations) { const option = create('option', configuration.label);option.value = configuration.id;option.dataset.url = new URL(configuration.url, url).href;configurations.append(option); }
    if ([...configurations.options].some(option => option.value === requested.get('configuration'))) configurations.value = requested.get('configuration');
    await load();
  })().catch(error => { root.replaceChildren(status);fail(error); });
})();
