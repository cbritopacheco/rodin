/* Browser regressions for a generated benchmark site; run against local build output. */
const assert = require('node:assert/strict');
const { pathToFileURL } = require('node:url');
const path = require('node:path');
const { chromium } = require('playwright');

(async () => {
  if (!process.argv[2]) throw new Error('Usage: node dev/test_benchmark_navigation.cjs SITE_DIRECTORY');
  const browser = await chromium.launch({
    headless: true,
    ...(process.env.RODIN_BROWSER_EXECUTABLE ? { executablePath: process.env.RODIN_BROWSER_EXECUTABLE } : {})
  });
  try {
    const firstLoad = await browser.newPage();
    await firstLoad.addInitScript(() => {
      window.domReady = false;
      document.addEventListener('DOMContentLoaded', () => { window.domReady = true; });
    });
    let releaseMath;
    await firstLoad.route('**/tex-svg.js', async route => {
      await new Promise(resolve => { releaseMath = resolve; });
      await route.fulfill({body: '', contentType: 'application/javascript'});
    });
    await firstLoad.goto(pathToFileURL(path.resolve(process.argv[2], 'index.html')).href, {waitUntil: 'commit'});
    await firstLoad.getByRole('button', {name: 'Color theme'}).waitFor();
    assert.equal(await firstLoad.evaluate(() => window.domReady), false);
    assert.ok(await firstLoad.locator('.rodin-theme-icon path').getAttribute('d'));
    await firstLoad.waitForFunction(() => document.querySelector('script[src*="tex-svg.js"]'));
    while (!releaseMath) await new Promise(resolve => setTimeout(resolve, 10));
    releaseMath();
    await firstLoad.waitForLoadState('domcontentloaded');
    assert.equal(await firstLoad.locator('#rodin-theme-trigger').count(), 1);
    await firstLoad.close();
    const page = await browser.newPage({ viewport: { width: 1280, height: 650 }, colorScheme: 'light' });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.goto(pathToFileURL(path.resolve(process.argv[2], 'index.html')).href);
    const navigation = page.getByRole('navigation', { name: 'On this page' });
    await navigation.waitFor();
    const content = page.locator('main > article > .m-container > .m-row > [data-rodin-content]');
    assert.ok(await content.evaluate(node => node.getBoundingClientRect().left > 200));
    const docs = await browser.newPage({ viewport: { width: 1280, height: 650 } });
    await docs.goto(pathToFileURL(path.resolve(process.argv[2], '../index.html')).href);
    const metrics = target => target.locator('h1').first().evaluate(node => ({
      x: node.getBoundingClientRect().left,
      width: node.getBoundingClientRect().width,
      size: getComputedStyle(node).fontSize
    }));
    const expected = await metrics(docs);
    const actual = await metrics(page);
    assert.ok(Math.abs(actual.x - expected.x) < 1);
    assert.ok(Math.abs(actual.width - expected.width) < 1);
    assert.equal(actual.size, expected.size);
    const searchCenter = (target, selector) => target.locator(selector).evaluate(node => {
      const bounds = node.getBoundingClientRect();
      return { x: bounds.left + bounds.width / 2, y: bounds.top + bounds.height / 2 };
    });
    const docsSearch = await searchCenter(docs, '#m-navbar-collapse .m-doc-search-icon svg');
    const benchmarkSearch = await searchCenter(page, '.benchmark-magnifier-icon');
    assert.ok(Math.abs(docsSearch.x - benchmarkSearch.x) < 1);
    assert.ok(Math.abs(docsSearch.y - benchmarkSearch.y) < 1);
    const back = await page.locator('#m-navbar-brand').boundingBox();
    const brand = await page.locator('.benchmark-search-trigger').boundingBox();
    assert.ok(back.x + back.width <= brand.x + 1);
    await docs.close();
    await navigation.hover();
    const control = page.getByRole('button', { name: 'Sidebar on right' });
    await control.click();
    assert.equal(await page.locator('body').getAttribute('data-contents-side'), 'right');
    assert.equal(await control.getAttribute('aria-pressed'), 'true');
    await navigation.hover();
    await page.getByRole('button', { name: 'Sidebar on left' }).click();
    assert.equal(await page.locator('body').getAttribute('data-contents-side'), 'left');
    const link = navigation.getByRole('link', { name: 'Infrastructure', exact: true });
    await link.click();
    await page.waitForFunction(() => {
      const heading = [...document.querySelectorAll('main h2')].find(node => node.textContent === 'Infrastructure');
      return heading && Math.abs(heading.getBoundingClientRect().top) < 145;
    });
    await page.waitForFunction(() => document.querySelector('.rodin-page-nav [aria-current="location"]')?.textContent === 'Infrastructure');
    // Manual scrolling, without another navigation click, must update the active section.
    await page.locator('main h3').filter({ hasText: 'Context' }).evaluate(node => window.scrollTo(0, window.scrollY + node.getBoundingClientRect().top - 120));
    await page.waitForFunction(() => document.querySelector('.rodin-page-nav [aria-current="location"]')?.textContent === 'Context');
    const focus = await page.evaluate(() => document.activeElement?.textContent);
    await page.mouse.wheel(0, 200);
    assert.equal(await page.evaluate(() => document.activeElement?.textContent), focus);
    await page.getByRole('button', { name: 'Color theme', exact: true }).click();
    await page.getByRole('menuitemradio', { name: 'Dark', exact: true }).click();
    assert.equal(await page.locator('html').getAttribute('data-theme'), 'dark');
    await page.goto(pathToFileURL(path.resolve(process.argv[2], 'Stokes/InfSupStability/index.html')).href);
    const spectrum = page.getByRole('navigation', { name: 'On this page' }).getByRole('link', { name: 'Generalized Schur spectrum', exact: true });
    await spectrum.click();
    await page.waitForFunction(() => document.querySelector('.rodin-page-nav [aria-current="location"]')?.textContent === 'Generalized Schur spectrum');
    await page.evaluate(() => window.MathJax?.startup?.promise);
    await page.setViewportSize({ width: 390, height: 800 });
    assert.ok(await page.locator('body').evaluate(node => node.scrollWidth <= window.innerWidth + 1));
    assert.deepEqual(errors, []);
    console.log('Benchmark navigation, destination icons, scrolling highlights and themes passed.');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
