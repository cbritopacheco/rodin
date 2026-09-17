# Custom agents for Rodin

This directory holds custom agent definitions for the Rodin finite element
framework. Each agent is a Markdown file with YAML frontmatter:

```markdown
---
name: <AgentName>
description: <When this agent should be selected>
---

<the agent's operating instructions>
```

## Available agents

### `rodin.agent.md` — Rodin

The contributor agent. It carries the repository's architecture rules, the
house style and the four CI style gates, how to scope builds and tests, the
verification standard, and the commit/PR conventions used here.

Use it for any code, test, or documentation change in this repository:

- `@Rodin add a Neumann term to the Poisson example and test it`
- `@Rodin why does the ClangTidy gate fail on my branch?`
- `@Rodin scope and run the tests affected by my Geometry change`

## How this relates to the other instruction files

Agent guidance is deliberately layered; none of these files should duplicate
another:

| File | Audience | Contains |
|---|---|---|
| `.github/agents/*.agent.md` | A selected custom agent | How to *work*: gates, scoping, verification, commit voice |
| `.github/copilot-instructions.md` | All Copilot requests in the repo | Architecture rules and style conventions in depth |
| `AGENTS.md` (repo root) | Every AI agent (Claude, Codex, Cursor, …) | Entry point, non-negotiable rules, pointers |
| `doc/agents/` | Every AI agent, on demand | The hierarchical knowledge base (philosophy, per-module, theory) |
| `CONTRIBUTING.md` | Humans and agents | The five style layers and how to run each check |

When a fact changes (a module moves, a gate is added, a command changes), update
the one file that owns it and let the others keep pointing.

## Environment

`.github/workflows/copilot-setup-steps.yml` provisions the toolchain and
dependencies for the GitHub coding agent. Keep it in step with the Build and
Tests workflows: an agent that cannot configure the project cannot verify its
own work.

## Adding an agent

1. Create `<name>.agent.md` here with `name` and `description` frontmatter.
2. Write instructions that are *specific to this repository* — general C++
   advice is noise. Prefer pointing at `doc/agents/` over restating it.
3. State the verification the agent owes: which tests, which gates, what
   evidence it must show.
4. Add it to the table above.
