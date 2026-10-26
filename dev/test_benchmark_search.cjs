/* Search uses the real pinned m.css engine and the generated benchmark index. */
const assert = require('node:assert/strict');
const { pathToFileURL } = require('node:url');
const path = require('node:path');
const { chromium } = require('playwright');

(async () => {
  if (!process.argv[2]) throw new Error('Usage: node dev/test_benchmark_search.cjs SITE_DIRECTORY');
  const browser = await chromium.launch({ headless: true,
    ...(process.env.RODIN_BROWSER_EXECUTABLE ? { executablePath: process.env.RODIN_BROWSER_EXECUTABLE } : {}) });
  try {
    const page = await browser.newPage({ viewport: { width: 1280, height: 850 }, colorScheme: 'light' });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.route('https://cdn.jsdelivr.net/npm/chart.js@*/**', route => route.fulfill({ contentType: 'application/javascript', body: '' }));
    await page.route('https://cdn.jsdelivr.net/npm/mathjax@*/**', route => route.fulfill({ contentType: 'application/javascript', body: '' }));
    const root = path.resolve(process.argv[2]);
    await page.goto(pathToFileURL(path.join(root, 'index.html')).href);
    await page.waitForFunction(() => window.Search?.trie);
    await page.getByRole('link', { name: 'Search benchmarks', exact: true }).click();
    const input = page.getByRole('searchbox', { name: 'Find benchmarks' });
    await input.fill('SpatialMatrix');
    const first = page.locator('#search-results > li').first();
    assert.match(await first.innerText(), /Modules\/Math\/SpatialMatrix/);
    assert.ok(await first.locator('.benchmark-search-brief').innerText());
    const urls = await page.locator('#search-results a').evaluateAll(links => links.map(link => link.href));
    assert.equal(urls.length, new Set(urls).size);
    await input.fill('Stokes stability');
    await page.locator('#search-results').waitFor({ state: 'visible' });
    assert.ok(await page.locator('#search-results a').count());
    await input.fill('zzzz-absent-benchmark');
    await page.locator('#search-notfound').waitFor({ state: 'visible' });
    await page.keyboard.press('Escape');
    await page.waitForFunction(() => location.hash !== '#search');
    assert.equal(await page.getByRole('link', { name: 'Search benchmarks', exact: true }).evaluate(node => node === document.activeElement), true);
    await page.goto(pathToFileURL(path.join(root, 'Stokes/InfSupStability/index.html')).href);
    await page.keyboard.press('t');
    await input.fill('Math::SpatialMatrix');
    assert.ok(await page.locator('#search-results a').count());
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowUp');
    await page.keyboard.press('Enter');
    await page.waitForURL('**/Rodin/Math/SpatialMatrix/index.html');
    await page.getByRole('link', { name: 'Search benchmarks', exact: true }).click();
    await input.fill('pyramid');
    await page.setViewportSize({ width: 390, height: 800 });
    assert.ok(await page.locator('#search').evaluate(node => node.scrollWidth <= window.innerWidth + 1));
    await page.getByRole('button', { name: 'Close benchmark search', exact: true }).click();
    assert.deepEqual(errors, []);
    console.log('Benchmark search: names, descriptions, scope aliases, deduplication, keyboard navigation, nested links and mobile layout passed.');
  } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
