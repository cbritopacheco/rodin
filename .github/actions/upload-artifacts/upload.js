'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');

function readInputs(env) {
  const number = (name, fallback, maximum) => {
    const value = env[name] || String(fallback);
    if (!/^\d+$/.test(value) || !Number.isSafeInteger(Number(value)) || Number(value) > maximum) {
      throw new Error(`Invalid upload option: ${name}`);
    }
    return Number(value);
  };
  const boolean = name => {
    const value = env[name] || 'false';
    if (!['true', 'false'].includes(value)) throw new Error(`Invalid upload option: ${name}`);
    return value === 'true';
  };
  const input = {
    name: env.UPLOAD_NAME, patterns: env.UPLOAD_PATH,
    overwrite: boolean('UPLOAD_OVERWRITE'), includeHidden: boolean('UPLOAD_HIDDEN'),
    compressionLevel: number('UPLOAD_COMPRESSION', 6, 9),
    retentionDays: number('UPLOAD_RETENTION', 0, Number.MAX_SAFE_INTEGER),
    noFiles: env.UPLOAD_NO_FILES || 'warn'
  };
  if (!input.name || !input.patterns || !['error', 'warn', 'ignore'].includes(input.noFiles)) {
    throw new Error('Expected artifact name, paths and a valid missing-file policy');
  }
  return input;
}

function transient(error) {
  const seen = new Set();
  let retry = false;
  for (let current = error; current && !seen.has(current); current = current.cause) {
    seen.add(current);
    const message = String(current.message || '');
    const status = Number(current.statusCode || current.status || current.response?.status ||
      message.match(/Failed request:\s*\((\d{3})\)/)?.[1]);
    if ([400, 401, 403, 404, 409, 422].includes(status) ||
        /quota|storage limit|permission denied|ENOSPC|EACCES|EINTEGRITY|ABORT_ERR/i.test(message) ||
        ['UsageError', 'AbortError', 'ArtifactNotFoundError', 'FilesNotFoundError'].includes(current.name)) return false;
    if ([408, 429, 500, 502, 503, 504].includes(status) ||
        /\b(?:ECONNRESET|ECONNREFUSED|ECONNABORTED|ETIMEDOUT|ESOCKETTIMEDOUT|ERR_SOCKET_TIMEOUT|EPIPE|ENOTFOUND|EAI_AGAIN|ENETUNREACH|EHOSTUNREACH)\b/.test(String(current.code) + ' ' + message) ||
        /\b(?:HTTP|status(?: code)?)[:\s]+(?:408|429|500|502|503|504)\b/i.test(message) ||
        message === 'Upload progress stalled.') retry = true;
  }
  return retry;
}

async function uploadArtifact(client, glob, input, core,
  sleep = milliseconds => new Promise(resolve => setTimeout(resolve, milliseconds))) {
  // Match upload-artifact's glob, hidden-file and least-common-root semantics.
  const globber = await glob.create(input.patterns, {followSymbolicLinks: true,
    implicitDescendants: true, omitBrokenSymbolicLinks: true, excludeHiddenFiles: !input.includeHidden});
  const files = [];
  for (const file of await globber.glob()) {
    if (!(await fs.stat(file)).isDirectory()) files.push(file);
  }
  if (!files.length) {
    if (input.noFiles === 'error') throw new Error('Required artifact has no matching files');
    core[input.noFiles === 'warn' ? 'warning' : 'info']('No matching files; optional artifact was not uploaded.');
    return;
  }
  const roots = globber.getSearchPaths();
  let root = roots[0];
  if (!root) throw new Error('Artifact glob has no root directory');
  for (const candidate of roots.slice(1)) {
    while (path.relative(root, candidate).split(path.sep)[0] === '..') {
      const parent = path.dirname(root);
      if (parent === root) throw new Error('Artifact paths have no common root');
      root = parent;
    }
  }
  if (files.length === 1 && roots.length === 1 && root === files[0]) root = path.dirname(root);
  const options = {compressionLevel: input.compressionLevel};
  if (input.retentionDays) options.retentionDays = input.retentionDays;
  const attempts = 3;
  for (let attempt = 1; attempt <= attempts; ++attempt) {
    try {
      if (input.overwrite) {
        try { await client.deleteArtifact(input.name); }
        catch (error) { if (error.name !== 'ArtifactNotFoundError') throw error; }
      }
      core.info(`Uploading ${input.name}: ${files.length} files, attempt ${attempt}/${attempts}`);
      const result = await client.uploadArtifact(input.name, files, root, options);
      if (!Number.isSafeInteger(result?.id) || result.id <= 0 ||
          !/^[a-f0-9]{64}$/i.test(result.digest || '')) {
        throw new Error('Artifact upload did not return a finalized ID and digest');
      }
      core.setOutput('artifact-id', result.id);
      core.setOutput('artifact-digest', result.digest);
      core.setOutput('artifact-url', `${process.env.GITHUB_SERVER_URL}/${process.env.GITHUB_REPOSITORY}/actions/runs/${process.env.GITHUB_RUN_ID}/artifacts/${result.id}`);
      return result;
    } catch (error) {
      if (attempt === attempts || !transient(error)) throw error;
      // No endpoint, token or signed blob URL belongs in retry messages.
      core.warning(`Transient artifact transfer failure for ${input.name}; retrying (${attempt + 1}/${attempts}).`);
      await sleep(2000 * attempt);
    }
  }
}

module.exports = {readInputs, transient, uploadArtifact};
