---
name: hvantk-workflow-conventions
description: The contract every hvantk workflow skill follows. Read first.
---

# hvantk workflow-skill conventions

A **workflow skill** tells a coding agent how to run hvantk for a user: which commands,
in what order, on which machine, with which checks. It never computes a result itself;
the hvantk commands it names are deterministic and the agent only orchestrates them.
A **resource spec** (`hvantk/skills/<provider>/SKILL.md`) is the other kind of skill in
this repository: it tells an agent how to *author* a data-source plugin. The two do not
overlap. A workflow skill may point at a resource spec; it never restates one.

## 1. Where things live

- `skills/<name>/SKILL.md`: the skill. Loaded whole into the agent's context on every
  turn after it is invoked, so it stays short.
- `skills/<name>/references/`: longer material, one file per step that needs it, read
  only when that step runs.
- `skills/<name>/templates/`: scripts the skill tells the agent to copy and run. Every
  site value in them comes from a variable, never a literal.
- `skills/<name>/eval/`: how the skill is measured (section 5). Removed from any sandbox
  in which an agent is evaluated.
- `skills/_site/`: the site profile (section 4).
- `skills/README.md`: the catalog. Its table and the folders must agree.

## 2. Frontmatter

Every `SKILL.md` opens with YAML frontmatter carrying exactly these keys:

| key | value |
| --- | --- |
| `name` | equals the folder name |
| `description` | starts with `Use when`; the phrases a harness matches an intent against |
| `stage` | `substrate`, `cohort`, `annotation` or `analysis` |
| `status` | `provisional` or `stable` |
| `scale` | `laptop`, `node` or `cluster`: the smallest setting the skill is meant for |
| `user-invocable` | `true` |
| `argument-hint` | one line |

## 3. The six sections

In this order, with these exact headings, each with a real body:

1. `## 1. Scope and outcome`
2. `## 2. Preflight`
3. `## 3. Procedure`
4. `## 4. Scale-up`
5. `## 5. Traps`
6. `## 6. Verification`

What goes where:

- **Scope and outcome.** What the workflow produces, when to use it, when not to.
- **Preflight.** Read the site profile; confirm the checkout commit and the container
  digest; confirm where inputs and outputs live; never guess a site fact. Stop and ask
  the user for any value that is missing.
- **Procedure.** Numbered steps. Each names the exact hvantk command or template it
  runs. Long material goes to `references/` and is named from the step.
- **Scale-up.** How the full-size run differs from a small one: resources, time,
  node-local scratch, checkpoints, fan-out.
- **Traps.** One line each: symptom, then cause, then fix.
- **Verification.** The commands a user runs on a real run, the values to expect, and
  the receipt the describe step writes, bound to the output's provenance sidecar.

Not in the body: the eval task, the reference values, any site value, anything
specific to one cohort or one publication. The checker rejects a body that names the
reference-values file.

## 4. The site profile

Every cluster-specific value comes from one shell-sourceable file that the user fills
in once and keeps out of git: copy `skills/_site/site-profile.example.sh` to
`skills/_site/site-profile.local.sh` and export `HVANTK_SITE_PROFILE` pointing at it.
`skills/_site/SITE-PROFILE.md` says what each variable means and how to find its value
on a new cluster. Templates `source "$HVANTK_SITE_PROFILE"` and refuse to run with a
required value missing (`skills/_site/check-site-profile.sh`). A skill that finds a
value missing stops and asks; it never guesses a partition, a path or an image.

## 5. Measured: the `eval/` folder

Five files, always:

| file | role |
| --- | --- |
| `task.md` | what an agent under evaluation is handed, verbatim. Names inputs, outputs, the budget, the job-name prefix and the `run-notes.md` it must write. Never names the reference values. |
| `stage.sh` | fetches and slices the public inputs once per cluster into `$HVANTK_EVAL_DATA/<name>/`, writing `staged.sha256`. Login-node work. |
| `control.sh` | the positive control: the skill's own procedure with no agent, then `check.py`, which must pass; and the negative control, a deliberately wrong output, on which `check.py --expect-fail` must pass. With `WRITE_REFERENCE=1` it (re)writes `reference.json`. |
| `check.py` | deterministic. Reads the output, `run-notes.md`, Slurm accounting and `reference.json`; writes `verdict.json` and `receipt.json`; exit 0 on pass. |
| `reference.json` | the values the control produced, with the commit, image digest and input digest that produced them. Never hand-edited. |

`run-notes.md` is three kinds of line, parsed by regular expression:

```
commit: <git sha, 7 to 40 hex>
container: <sha256 of the image file>
job: <slurm job id>        (one line per job submitted)
```

`verdict.json` carries five graded stages and the pass flag:

| stage | passes when |
| --- | --- |
| M1 | `run-notes.md` present; commit and container digest equal the expected ones |
| M2 | at least one Slurm job carries the run's name prefix; its partition is one the site profile names; CPUs, memory and time limit are within the task budget |
| M3 | such a job reached `COMPLETED`, and every job id in the notes exists in accounting |
| M4 | the output exists and passes its dataset contract |
| M5 | the output's values match `reference.json` within the stated tolerance |

`pass` is M4 and M5. A task that needs no job records M2 and M3 as `null`.

## 6. Sensitive data

No usernames, home or work paths, hostnames, partition, account or module names, e-mail
addresses or job ids anywhere under `skills/`. `hvantk/tests/test_no_site_data_in_skills.py`
rejects the generic shapes; the author's publish guard rejects the specific ones. Write
`$HVANTK_WORK`, `$HVANTK_CPU_PARTITION`, `<cluster-user>`, `you@example.org`.

## 7. Template

Copy this block verbatim into `skills/<name>/SKILL.md` and replace the `<...>`
placeholders. The headings are matched by whole-line equality; do not annotate them. Set
`stage` and `scale` to your skill's values; the template shows `substrate` and `node`.
`hvantk/tests/test_workflow_skill_conformance.py` checks the template itself.

```markdown
---
name: <folder-name>
description: Use when <the intents this skill serves, in the user's words>.
stage: substrate
status: provisional
scale: node
user-invocable: true
argument-hint: "<what to pass>"
---

# <Title>

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

<What this produces. When to use it. When not to.>

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE"` and run `bash skills/_site/check-site-profile.sh`; stop and ask for any missing value.
2. Checkout: confirm the commit the user expects is checked out on the machine that runs the job.
3. Image: confirm the image path and record its digest.
4. Inputs and outputs: confirm where they live and that the output location is writable from compute nodes.

## 3. Procedure

1. <Step, naming the exact command or template.>
2. <Step.>

## 4. Scale-up

<How the full-size run differs: resources, time, scratch, fan-out.>

## 5. Traps

- <Symptom.> Cause: <cause>. Fix: <fix>.

## 6. Verification

<The commands to run on a real run, the values to expect, and the receipt.>
```
