/* Test the published-history UI contract; numerical data and Chart rendering are fixtures. */
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const http = require('node:http');
const path = require('node:path');
const { chromium } = require('playwright');

(async () => {
  if (!process.argv[2]) throw new Error('Usage: node dev/test_benchmark_charts.cjs SITE_DIRECTORY');
  const root = path.resolve(process.argv[2]);
  const caseHistory = (name, metric, value, unit) => ({ schema: 1, window: 20, format: 'googlecpp', name, metric, entries: [
    { commit: 'a'.repeat(40), date: '2026-01-01', run: '1.1', results: [{ name, metric, value, unit, samples: [value] }] }
  ] });
  const histories = {
    '/__fixture/index.json': { schema: 1, configurations: [
      { id: 'cpu', label: 'CPU fixture', url: 'cpu.json' },
      { id: 'wall', label: 'Wall fixture', url: 'wall.json' }
    ] },
    '/__fixture/cpu.json': { schema: 1, cases: [
      { id: 'three', name: 'Kernel/3', metric: 'cpu_time', url: 'three.json' },
      { id: 'four', name: 'Kernel/4', metric: 'cpu_time', url: 'four.json' }
    ] },
    '/__fixture/wall.json': { schema: 1, cases: [
      { id: 'wall-three', name: 'Kernel/3', metric: 'real_time', url: 'wall-three.json' }
    ] },
    '/__fixture/three.json': caseHistory('Kernel/3', 'cpu_time', 2, 'us'),
    '/__fixture/four.json': caseHistory('Kernel/4', 'cpu_time', 5, 'ns'),
    '/__fixture/wall-three.json': caseHistory('Kernel/3', 'real_time', 9, 'ns')
  };
  const server = http.createServer(async (request, response) => {
    try {
      const url = new URL(request.url, 'http://localhost').pathname;
      if (histories[url]) { response.setHeader('Content-Type', 'application/json'); response.end(JSON.stringify(histories[url])); return; }
      const file = path.resolve(root, '.' + decodeURIComponent(url));
      if (!file.startsWith(root + path.sep)) throw new Error('Outside fixture site');
      const extension = path.extname(file);
      response.setHeader('Content-Type', { '.html': 'text/html', '.js': 'application/javascript', '.css': 'text/css' }[extension] || 'application/octet-stream');
      let content = await fs.readFile(file);
      if (url.endsWith('/SpatialMatrix/index.html')) content = content.toString().replace('data-benchmark-history', 'data-benchmark-history data-history-index="/__fixture/index.json"');
      response.end(content);
    } catch (_error) { response.statusCode = 404; response.end(); }
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const browser = await chromium.launch({ headless: true,
    ...(process.env.RODIN_BROWSER_EXECUTABLE ? { executablePath: process.env.RODIN_BROWSER_EXECUTABLE } : {}) });
  try {
    const page = await browser.newPage();
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.route('https://cdn.jsdelivr.net/npm/chart.js@*/**', route => route.fulfill({ contentType: 'application/javascript', body:
      'window.Chart=class{constructor(canvas,config){window.lastChart=config;}destroy(){}};' }));
    await page.route('https://cdn.jsdelivr.net/npm/mathjax@*/**', route => route.fulfill({ contentType: 'application/javascript', body: '' }));
    await page.goto(`http://127.0.0.1:${server.address().port}/Rodin/Math/SpatialMatrix/index.html`);
    await page.waitForFunction(() => window.lastChart?.data.datasets[0].data[0] === 2000);
    const filter = page.getByRole('searchbox', { name: 'Filter benchmark cases' });
    await filter.fill('Kernel/4');
    await page.waitForFunction(() => document.querySelector('[aria-label="Benchmark case"]').options.length === 1);
    await filter.pressSequentially('t');
    await page.waitForFunction(() => document.querySelector('[aria-label="Benchmark case"]').disabled);
    assert.notEqual(await page.evaluate(() => location.hash), '#search');
    await filter.fill('Kernel/4');
    await page.waitForFunction(() => window.lastChart.data.datasets[0].data[0] === 5);
    await filter.fill('');
    await page.waitForFunction(() => document.querySelector('[aria-label="Benchmark case"]').options.length === 2);
    await page.getByLabel('Benchmark configuration', { exact: true }).selectOption('wall');
    await page.waitForFunction(() => window.lastChart.data.datasets[0].data[0] === 9);
    assert.match(page.url(), /configuration=wall/);
    assert.match(page.url(), /case=wall-three/);
    await page.reload();
    await page.waitForFunction(() => window.lastChart?.data.datasets[0].data[0] === 9);
    assert.equal(await page.getByLabel('Benchmark configuration', { exact: true }).inputValue(), 'wall');
    assert.match(await page.evaluate(() => window.lastChart.options.scales.yAxes[0].scaleLabel.labelString), /real_time/);
    await page.setViewportSize({ width: 1280, height: 900 });
    const configBox = await page.getByLabel('Benchmark configuration', { exact: true }).boundingBox();
    const filterBox = await filter.boundingBox();
    const caseBox = await page.getByLabel('Benchmark case', { exact: true }).boundingBox();
    assert.ok(filterBox.x >= configBox.x + configBox.width + 12);
    assert.ok(caseBox.y >= configBox.y + configBox.height + 12);
    assert.ok(caseBox.width > configBox.width * 1.8);
    assert.equal(await page.locator('.benchmark-selected-case').textContent(), 'Kernel/3');
    await page.setViewportSize({ width: 390, height: 800 });
    const mobileConfig = await page.getByLabel('Benchmark configuration', { exact: true }).boundingBox();
    const mobileFilter = await filter.boundingBox();
    assert.ok(mobileFilter.y >= mobileConfig.y + mobileConfig.height + 12);

    assert.ok(await page.locator('body').evaluate(node => node.scrollWidth <= window.innerWidth + 1));
    assert.deepEqual(errors, []);
    console.log('History loading, configuration/case selection, unit conversion and responsive controls passed.');
  } finally { await browser.close(); await new Promise(resolve => server.close(resolve)); }
})().catch(error => { console.error(error); process.exitCode = 1; });
