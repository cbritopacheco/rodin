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
    for (const option of document.querySelectorAll('.rodin-theme-menu [role="menuitemradio"]')) {
      option.setAttribute('aria-checked', String(option.dataset.theme === theme));
    }
    const light = theme === 'light' || (theme === 'system' && system.matches);
    document.querySelector('meta[name="theme-color"]')?.setAttribute('content', light ? '#f7f8fa' : '#20252d');
  }
  applyTheme();
  system.addEventListener('change', applyTheme);

  function addNavigationIcons() {
    const paths = {
      'pages.html': 'M6 3h9l4 4v14H6z M14 3v5h5 M9 12h7 M9 16h7',
      'modules.html': 'm12 3 9 5-9 5-9-5z M3 12l9 5 9-5 M3 16l9 5 9-5',
      'annotated.html': 'M4 3h16v7H4z M8 6h8 M12 10v4 M5 14h14 M5 14v3 M19 14v3 M2 17h6v4H2z M16 17h6v4h-6z',
      'files.html': 'M6 3h9l4 4v14H6z M14 3v5h5'
    };
    for (const link of document.querySelectorAll('#m-navbar-collapse > .m-row > ol > li > a')) {
      const drawing = paths[link.getAttribute('href')];
      if (!drawing || link.querySelector('svg')) continue;
      const icon = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      icon.setAttribute('viewBox', '0 0 24 24');
      icon.setAttribute('aria-hidden', 'true');
      icon.classList.add('rodin-nav-icon');
      const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      path.setAttribute('d', drawing);
      icon.append(path);
      link.prepend(icon);
    }
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', addNavigationIcons, {once: true});
  } else addNavigationIcons();

  function initializeControl() {
    const brand = document.querySelector('#m-navbar-brand');
    if (!brand) return false;
    if (document.querySelector('#rodin-theme-trigger')) return true;
    if (!document.body?.classList.contains('rodin-benchmarks')) {
      const home = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      home.setAttribute('viewBox', '2 2 20 20');
      home.setAttribute('aria-hidden', 'true');
      home.classList.add('rodin-home-icon');
      const homePath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      homePath.setAttribute('d', 'M3 10 12 3l9 7 M5 9v12h5v-7h4v7h5V9');
      home.append(homePath);
      brand.prepend(home);
    }
    const github = document.querySelector('#m-navbar-collapse a[href="https://github.com/cbritopacheco/rodin"]');
    if (github) {
      const mark = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      mark.setAttribute('viewBox', '0 0 16 16');
      mark.setAttribute('aria-hidden', 'true');
      mark.classList.add('rodin-github-icon');
      const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      path.setAttribute('d', 'M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82A7.65 7.65 0 0 1 8 3.86c.68 0 1.36.09 2 .26 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z');
      mark.append(path);
      github.prepend(mark);
    }
    const item = document.createElement('div');
    item.className = 'rodin-theme-control';
    const trigger = document.createElement('button');
    trigger.type = 'button';
    trigger.id = 'rodin-theme-trigger';
    trigger.className = 'rodin-theme-trigger';
    trigger.setAttribute('aria-label', 'Color theme');
    trigger.setAttribute('aria-haspopup', 'menu');
    trigger.setAttribute('aria-expanded', 'false');
    trigger.setAttribute('aria-controls', 'rodin-theme-menu');
    const icon = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    icon.setAttribute('viewBox', '0 0 24 24');
    icon.setAttribute('aria-hidden', 'true');
    icon.classList.add('rodin-theme-icon');
    icon.append(document.createElementNS('http://www.w3.org/2000/svg', 'path'));
    trigger.append(icon);
    const menu = document.createElement('div');
    menu.id = 'rodin-theme-menu';
    menu.className = 'rodin-theme-menu';
    menu.setAttribute('role', 'menu');
    menu.setAttribute('aria-labelledby', trigger.id);
    menu.hidden = true;
    const options = [];

    function closeMenu(focus = false) {
      menu.hidden = true;
      trigger.setAttribute('aria-expanded', 'false');
      if (focus) trigger.focus();
    }

    function openMenu(index = options.findIndex(option => option.dataset.theme === theme)) {
      menu.hidden = false;
      trigger.setAttribute('aria-expanded', 'true');
      options[index].focus();
    }

    for (const [value, text] of [['system', 'System'], ['light', 'Light'], ['dark', 'Dark']]) {
      const option = document.createElement('button');
      option.type = 'button';
      option.setAttribute('role', 'menuitemradio');
      option.tabIndex = -1;
      option.dataset.theme = value;
      option.textContent = text;
      option.addEventListener('click', () => {
        theme = value;
        try { localStorage.setItem(key, theme); } catch {}
        const url = new URL(location.href);
        url.searchParams.delete('theme');
        history.replaceState(history.state, '', url);
        applyTheme();
        closeMenu(true);
      });
      options.push(option);
      menu.append(option);
    }
    trigger.addEventListener('click', () => menu.hidden ? openMenu() : closeMenu());
    trigger.addEventListener('keydown', event => {
      if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
        event.preventDefault();
        openMenu(event.key === 'ArrowDown' ? 0 : options.length - 1);
      }
    });
    menu.addEventListener('keydown', event => {
      const index = options.indexOf(document.activeElement);
      if (event.key === 'Escape') {
        event.preventDefault();
        closeMenu(true);
      } else if (event.key === 'Tab') {
        closeMenu(true);
      } else if (['ArrowDown', 'ArrowUp', 'Home', 'End'].includes(event.key)) {
        event.preventDefault();
        const next = event.key === 'Home' ? 0 : event.key === 'End' ? options.length - 1
          : (index + (event.key === 'ArrowDown' ? 1 : options.length - 1)) % options.length;
        options[next].focus();
      }
    });
    document.addEventListener('pointerdown', event => {
      if (!item.contains(event.target)) closeMenu();
    });
    item.addEventListener('focusout', event => {
      if (!item.contains(event.relatedTarget)) closeMenu();
    });
    item.append(trigger, menu);
    brand.closest('nav').append(item);
    applyTheme();
    return true;
  }

  if (!initializeControl()) {
    const observer = new MutationObserver(() => {
      if (initializeControl()) observer.disconnect();
    });
    observer.observe(document.documentElement, {childList: true, subtree: true});
    document.addEventListener('DOMContentLoaded', () => {
      initializeControl();
      observer.disconnect();
    }, {once: true});
  }
})();
