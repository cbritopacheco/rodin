# GitHub Copilot Agents for Rodin

This directory contains custom GitHub Copilot agents that provide specialized assistance for the Rodin finite element framework.

## Available Agents

### Rodin Agent

**File**: `rodin.yml`

**Purpose**: The Rodin agent is an expert at building and testing Rodin code changes. It automatically compiles the codebase and runs tests whenever modifications are detected.

**Capabilities**:
- Compiles the Rodin codebase using CMake
- Runs unit tests, manufactured tests, and benchmarks
- Reports build and test results clearly
- Suggests fixes for compilation errors or test failures
- Provides guidance on the build system and test infrastructure

**When to Use**:
- After making changes to source code
- When you need to verify tests pass
- To troubleshoot build or test failures
- For guidance on the CMake build system
- To understand test results and failures

**How to Use**:
In GitHub Copilot Chat, you can invoke the Rodin agent with commands like:
- `@Rodin build and test my changes`
- `@Rodin compile the code and run unit tests`
- `@Rodin why did the tests fail?`
- `@Rodin help me fix this build error`

## Agent Structure

Each agent is defined in a YAML file with the following structure:

```yaml
name: <AgentName>
description: <Brief description of the agent's purpose>
instructions: |
  <Detailed instructions for the agent's behavior>
```

## Creating New Agents

To create a new custom agent:

1. Create a new YAML file in this directory (e.g., `my-agent.yml`)
2. Define the agent's name, description, and instructions
3. Follow the structure of existing agents
4. Document the agent in this README

## Best Practices

- Keep agent instructions clear and specific
- Include examples of commands and workflows
- Document common issues and troubleshooting steps
- Provide context about the Rodin codebase
- Test agent responses to ensure they work as expected

## References

- [Rodin Build System](../../CMakeLists.txt)
- [Copilot Instructions](../copilot-instructions.md)
- [Testing Documentation](../../tests/README.md)
- [GitHub Copilot Documentation](https://docs.github.com/en/copilot)
