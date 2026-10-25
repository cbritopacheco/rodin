(() => {
  const key = 'rodin-color-theme';
  const system = window.matchMedia('(prefers-color-scheme: light)');
  const requested = new URLSearchParams(location.search).get('theme');
  let saved = null;
  try { saved = localStorage.getItem(key); } catch {}
  let theme = ['light', 'dark', 'system'].includes(requested) ? requested
    : ['light', 'dark', 'system'].includes(saved) ? saved : 'system';

  const icons = {
    system: 'M3 4h18v13H3z M8 21h8 M12 17v4',
    light: 'M16 12a4 4 0 1 1-8 0 4 4 0 0 1 8 0 M12 2v2 M12 20v2 M2 12h2 M20 12h2 M5 5l1.5 1.5 M17.5 17.5L19 19 M5 19l1.5-1.5 M17.5 6.5L19 5',
    dark: 'M21 13a9 9 0 0 1-10-10A9 9 0 1 0 21 13Z'
  };

  function applyTheme() {
    if (theme === 'system') delete document.documentElement.dataset.theme;
    else document.documentElement.dataset.theme = theme;
    document.querySelector('.rodin-theme-icon path')?.setAttribute('d', icons[theme]);
    document.querySelector('[aria-label="Color theme"]')?.setAttribute('title', 'Color theme: ' + theme);
    const light = theme === 'light' || (theme === 'system' && system.matches);
    document.querySelector('meta[name="theme-color"]')?.setAttribute('content', light ? '#f7f8fa' : '#20252d');
  }
  applyTheme();
  system.addEventListener('change', applyTheme);

  document.addEventListener('DOMContentLoaded', () => {
    const brand = document.querySelector('#m-navbar-brand');
    if (!brand) return;
    const item = document.createElement('div');
    item.className = 'rodin-theme-control';
    const label = document.createElement('label');
    label.title = 'Choose color theme';
    const icon = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    icon.setAttribute('viewBox', '0 0 24 24');
    icon.setAttribute('aria-hidden', 'true');
    icon.classList.add('rodin-theme-icon');
    icon.append(document.createElementNS('http://www.w3.org/2000/svg', 'path'));
    label.append(icon);
    const select = document.createElement('select');
    select.setAttribute('aria-label', 'Color theme');
    for (const [value, text] of [['system', 'System'], ['light', 'Light'], ['dark', 'Dark']]) {
      const option = document.createElement('option');
      option.value = value;
      option.textContent = text;
      select.append(option);
    }
    select.value = theme;
    select.addEventListener('change', () => {
      theme = select.value;
      try { localStorage.setItem(key, theme); } catch {}
      const url = new URL(location.href);
      url.searchParams.delete('theme');
      history.replaceState(history.state, '', url);
      applyTheme();
    });
    label.append(select);
    item.append(label);
    brand.closest('nav').append(item);
    applyTheme();
  });
})();
