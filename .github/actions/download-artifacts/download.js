'use strict';

const fs = require('node:fs/promises');
const path = require('node:path');

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

async function downloadArtifacts(client, names, destination, core,
  sleep = milliseconds => new Promise(resolve => setTimeout(resolve, milliseconds))) {
  names = await names;
  await fs.mkdir(destination, {recursive: true});
  // Exact runtime-service lookups bypass both the 1000-item listing cap and
  // the installation REST API quota. Never supply findBy (a REST API request).
  for (const name of names) {
    const target = path.join(destination, name);
    const staging = await fs.mkdtemp(path.join(destination, '.download-'));
    try {
      const {artifact} = await client.getArtifact(name);
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
      await fs.rm(target, {recursive: true, force: true});
      await fs.rename(staging, target);
      core.info(`Published mandatory artifact ${name} to ${target}`);
    } finally {
      await fs.rm(staging, {recursive: true, force: true});
    }
  }
}

module.exports = {readNames, downloadArtifacts};
