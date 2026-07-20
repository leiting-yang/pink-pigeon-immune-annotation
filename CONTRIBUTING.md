# Contributing to the Workflow

First off, thank you for considering contributing to this project! Your contributions help make it better for everyone in the lab.

## Guiding Philosophy

This repository exists to **standardize workflows** so that anyone in the lab—from master's students to postdocs—can run analyses within days, not weeks.

We recognize that:
- **Bug fixes often emerge from real research needs**: You're using the code for your own project, you hit a snag, and you adapt it to make it work for your case.
- **Not everyone has the same software development background**: Some fixes are "good enough for me" rather than "production-ready." That's okay.

**Our goal is to move from "it works for me" to "it works for everyone."** The GitHub repository is where we collect and refine those improvements collaboratively.

## Two Ways to Contribute

Because different people work in different ways, we support **two strategies**. Choose whichever feels most natural to you:

### Option A: Fork the Repository (For Personal Experimentation)

If you prefer to experiment freely before contributing:

1. Fork this repository to your personal GitHub account.
2. Clone your fork and create a branch:
```bash
git clone https://github.com/your-username/immune-annotation-workflow.git
git checkout -b feature/your-feature-name
```
1. Make your changes, commit, and push to your fork.
2. Open a Pull Request from your branch to the `main` branch of the organization repository.

**Important**: If your fix is purely for your own use case and not intended for the broader lab, you can keep it in your fork. But if you think it could benefit others, **please open a PR** so we can discuss and integrate it.

### Option B: Work Directly on the Organization Repository (Recommended for Active Contributors)

If you have write access to this repository, you can:

1. Clone the repository directly:
```bash
git clone https://github.com/EvoConsGen/workflow.git
```
2. Create a branch for your changes:
```bash
git checkout -b feature/your-feature-name
```
3. Make your changes, commit them, and push:
```bash
git push origin feature/your-feature-name
```
4. Open a Pull Request (PR) against the `main` branch. 

## Branch Naming Conventions (Suggested, Not Mandatory)

While we don't enforce strict naming rules, using descriptive branch names helps everyone understand what's being worked on. Please consider using:

- `bugfix/short-description` for bug fixes
- `feature/short-description` for new functionality
- `docs/short-description` for documentation updates
- `experiment/short-description` for exploratory changes

Example: `bugfix/fix-input-parsing`, `feature/add-faster-aligner`

## Commit Messages

Write clear, concise commit messages that explain **what** changed and **why**. A good format is:

```
Short summary (50 chars or less)

Optional longer explanation with context.
```

Example:
```
Add support for new aligner v2.0

The new aligner is 30% faster and produces better results for
long reads. Updated the config file and added a flag.
```

## Pull Request Process

We use **Pull Requests (PRs)** as the primary mechanism for integrating changes into the `main` branch. This ensures code quality and shared knowledge.

### Step 1: Open a PR
- Go to the repository and click "New Pull Request."
- Select your branch and compare it to `main`.
- Write a descriptive title and explain what your changes do and why they're needed.

### Step 2: Review and Approval
- **One approval is required** before merging.
- The **EvoConsGen Review Team** (a designated GitHub team) will review your PR. Other team members may also review.
- All comments and conversations must be **resolved** before merging. If there's a discussion, we expect it to reach a conclusion.

### Step 3: Merging
- **Once your PR has approval and all conversations are resolved, you (the author) are responsible for merging it.**
- You can use **"Merge"** or **"Squash and merge"** (we recommend squash for cleaner history).
- **Do NOT use "Rebase and merge"** or `git push --force` on shared branches. This can rewrite history and confuse others.

## Code Style

We don't enforce a specific linter or formatter, but we do ask that you **maintain stylistic consistency** with the existing code. If you're adding new functions or modifying existing ones:

- Follow the naming conventions already used (e.g., snake_case for Python, camelCase for R if that's the style).
- Keep the code readable and add comments where logic is non-obvious.
- If you're unsure, mimic the style of the most recent code in the repository.

## Documentation

Good documentation is essential for a workflow that others will use.

- **For new features**: If you add a new functionality, please update the README with a section explaining how to use it.
- **For bug fixes**: If the fix changes how something works, update the relevant docstrings or comments.
- **Before implementing a major new feature**: Open an Issue or start a Discussion to get feedback from the team. This helps avoid duplicated effort and ensures the feature aligns with the lab's needs.

## Reporting Bugs and Suggesting Improvements

We use **GitHub Issues** for tracking bugs and feature requests.

- **No strict template** — write in free text, but please include:
	- What you were trying to do
	- What actually happened (error messages, unexpected outputs)
	- Steps to reproduce the issue (if possible)
	- Your environment (OS, software versions)

- If you're unsure whether something is a bug or a feature, open a **Discussion** first and we can move it to an Issue if needed.

## Communication and Collaboration

- For **general questions** about the workflow, use the **Discussions** tab.
- For **specific technical issues**, open an **Issue**.
- For **proposing changes**, open a **PR** (or discuss via Issue first for larger changes).

We encourage everyone to participate in reviews — reviewing others' code is a great way to learn and share knowledge across the lab.

## Acknowledgment

Contributions are recognized in the repository's `AUTHORS.md` file. If you contribute significantly, feel free to add your name or ask the maintainer to do so.

## Questions?

If anything in this guide is unclear, or if you're stuck on the technical process, don't hesitate to reach out to the repository maintainer or open a Discussion.

**Thank you for contributing!** Every fix, every improvement, and every shared insight makes this workflow more useful for everyone in the lab.
