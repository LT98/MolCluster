# Archive

Finished work. **Nothing here is scheduled**, and nothing here should be read to find out what
to do next — read it to find out why something is the way it is.

| File | Holds | Active counterpart |
|---|---|---|
| [`BUGS_resolved.md`](BUGS_resolved.md) | defects that are fixed and covered by a test, with the diagnosis that produced the fix | [`../BUGS.md`](../BUGS.md) |
| [`PLAN_completed.md`](PLAN_completed.md) | milestones M0–M4 with their exit gates, the decided interface architecture, the `legacy/` salvage ledger, and the full changelog rev 1 → 23 | [`../PLAN_implementation.md`](../PLAN_implementation.md) |
| [`DESIGN_history.md`](DESIGN_history.md) | every revision of the design doc, newest first | [`../DESIGN_registry_assembly.md`](../DESIGN_registry_assembly.md) |

**The rule that keeps this useful:** an item moves here the moment it is done. An active
document never carries a struck-through entry or a section marked ✅ — its length is supposed to
be the size of the remaining problem.

**This is also where bug-fix rationale lives** rather than in source comments. If you fixed
something and the reasoning is worth keeping, it goes in `BUGS_resolved.md`; the code gets a
line naming the invariant, not the history.
