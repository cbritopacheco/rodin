'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');
const os = require('node:os');

async function readNames(environment) {
  const inline = environment.ARTIFACT_NAMES;
  const file = environment.ARTIFACT_NAMES_FILE;
  if (Boolean(inline) === Boolean(file)) {
    throw new Error('Specify exactly one artifact name list or list file');
  }
  const names = JSON.parse(inline || await fs.readFile(file, 'utf8'));
  if (!Array.isArray(names) || !names.length ||
      names.some(name => typeof name !== 'string' || !/^[A-Za-z0-9_-]+$/.test(name)) ||
      new Set(names).size !== names.length) {
    throw new Error('Expected a nonempty list of unique exact artifact names');
  }
  return names;
}

function transientLookup(error) {
  const seen = new Set();
  let retry = false;
  for (let current = error; current && !seen.has(current); current = current.cause) {
    seen.add(current);
    const message = String(current.message || '');
    const status = Number(current.statusCode || current.status || current.response?.status ||
      message.match(/Failed request:\s*\((\d{3})\)/)?.[1]);
    if ([400, 401, 403, 404, 409, 422].includes(status) ||
        ['ArtifactNotFoundError', 'UsageError', 'AbortError'].includes(current.name) ||
        /quota|permission denied|not found|ABORT_ERR/i.test(message)) return false;
    if ([408, 429, 500, 502, 503, 504].includes(status) ||
        /\b(?:ECONNRESET|ECONNREFUSED|ECONNABORTED|ETIMEDOUT|ESOCKETTIMEDOUT|EPIPE|ENOTFOUND|EAI_AGAIN|ENETUNREACH|EHOSTUNREACH)\b/.test(String(current.code) + ' ' + message) ||
        /\b(?:HTTP|status(?: code)?)[:\s]+(?:408|429|500|502|503|504)\b/i.test(message)) retry = true;
  }
  return retry;
}

async function downloadArtifacts(client, names, destination, core,
  sleep = milliseconds => new Promise(resolve => setTimeout(resolve, milliseconds)), mergeSingle = false) {
  names = await names;
  if (mergeSingle && names.length !== 1) {
    throw new Error('Merged extraction requires exactly one artifact');
  }
  await fs.mkdir(destination, {recursive: true});
  // Exact runtime-service lookups bypass both the 1000-item listing cap and
  // the installation REST API quota. Never supply findBy (a REST API request).
  for (const name of names) {
    const target = path.join(destination, name);
    const staging = await fs.mkdtemp(path.join(mergeSingle ? os.tmpdir() : destination, '.download-'));
    try {
      let artifact;
      for (let attempt = 1; ; ++attempt) {
        try {
          ({artifact} = await client.getArtifact(name));
          break;
        } catch (error) {
          if (attempt === 3 || !transientLookup(error)) throw error;
          core.warning(`Transient artifact lookup failure for ${name}; retrying (${attempt + 1}/3).`);
          await sleep(2000 * attempt);
        }
      }
      if (artifact.name !== name) {
        throw new Error(`Artifact lookup returned the wrong name for ${name}`);
      }
      const attempts = 3;
      for (let attempt = 1; ; ++attempt) {
        try {
          core.info(`Downloading ${name} (ID ${artifact.id}, attempt ${attempt})`);
          const result = await client.downloadArtifact(artifact.id, {
            path: staging, ...(artifact.digest ? {expectedHash: artifact.digest} : {})
          });
          if (result?.digestMismatch) {
            throw new Error(`Artifact digest mismatch: ${name}`);
          }
          if (!(await fs.readdir(staging)).length) {
            throw new Error(`Empty artifact download: ${name}`);
          }
          break;
        } catch (error) {
          if (attempt === attempts) throw error;
          core.warning(`Retrying artifact ID ${artifact.id}: ${error.message}`);
          await fs.rm(staging, {recursive: true, force: true});
          await fs.mkdir(staging);
          await sleep(2000 * attempt);
        }
      }
      if (mergeSingle) {
        // Existing workspace files survive. Only fully verified contents merge.
        await fs.cp(staging, destination, {recursive: true});
      } else {
        await fs.rm(target, {recursive: true, force: true});
        await fs.rename(staging, target);
      }
      core.info(`Published mandatory artifact ${name} to ${mergeSingle ? destination : target}`);
    } finally {
      await fs.rm(staging, {recursive: true, force: true});
    }
  }
}

module.exports = {readNames, downloadArtifacts};
