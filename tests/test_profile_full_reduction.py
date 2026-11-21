"""Profiler measurements preserve calls and do not double-count nested work."""
from types import SimpleNamespace
import cProfile

import pytest

from scripts import profile_full_reduction as profiling


@pytest.mark.parametrize('same_category', [False, True])
def test_nested_stage_seconds_compose_without_overlap(monkeypatch, same_category):
    ticks = iter([0., 1., 3., 10.])
    monkeypatch.setattr(profiling, 'perf_counter', lambda: next(ticks))
    stages = profiling.StageTimes()
    token = object()
    owner = SimpleNamespace(child=lambda value: value)
    owner.parent = lambda value: owner.child(value)
    stages.wrap(owner, 'child', 'aperture' if same_category else 'psf')
    stages.wrap(owner, 'parent', 'aperture')
    assert owner.parent(token) is token
    assert sum(row['exclusive_seconds'] for row in stages.groups.values()) == 10.
    assert stages.groups['aperture']['inclusive_seconds'] == 10.
    if same_category:
        assert stages.groups['aperture']['calls'] == 2
    else:
        assert stages.groups['psf']['exclusive_seconds'] == 2.
        assert stages.groups['aperture']['exclusive_seconds'] == 8.
    assert stages.stack == []


def test_stage_failures_propagate_and_leave_usable_measurements(monkeypatch):
    ticks = iter([0., 2.])
    monkeypatch.setattr(profiling, 'perf_counter', lambda: next(ticks))
    stages = profiling.StageTimes()
    def fail():
        raise RuntimeError('actual reduction failure')
    owner = SimpleNamespace(operation=fail)
    stages.wrap(owner, 'operation', 'aperture')
    with pytest.raises(RuntimeError, match='actual reduction failure'):
        owner.operation()
    assert stages.groups['aperture']['exclusive_seconds'] == 2.
    assert stages.stack == []


def test_raw_call_export_preserves_import_aliases_with_identical_source_labels():
    first, second = {}, {}
    code = 'def evaluate():\n    return 1\n'
    exec(compile(code, 'same_source.py', 'exec'), first)
    exec(compile(code, 'same_source.py', 'exec'), second)
    profiler = cProfile.Profile()
    def caller():
        first['evaluate']()
        second['evaluate']()
        second['evaluate']()
    profiler.enable()
    caller()
    profiler.disable()
    rows = [row for row in profiling.raw_profile_rows(profiler) if row['file'] == 'same_source.py']
    assert sorted(row['calls'] for row in rows) == [1, 2]
    assert len({row['code_id'] for row in rows}) == 2
    all_rows = profiling.raw_profile_rows(profiler)
    ids = {row['code_id'] for row in all_rows}
    assert any(row['callees'] for row in all_rows)
    assert all(child['code_id'] in ids for row in all_rows for child in row['callees'])


def test_raw_call_export_matches_recreated_builtin_name_strings():
    name='<built-in method example>'
    edge_name=name.encode().decode()
    assert edge_name == name and edge_name is not name
    def entry(code,children=None):
        return SimpleNamespace(code=code,callcount=1,reccallcount=0,
                               inlinetime=.1,totaltime=.2,calls=children)
    profiler=SimpleNamespace(getstats=lambda:[entry(name),entry(test_raw_call_export_matches_recreated_builtin_name_strings.__code__,[entry(edge_name)])])
    rows=profiling.raw_profile_rows(profiler)
    assert rows[1]['callees'][0]['code_id']==rows[0]['code_id']


def test_raw_call_export_retains_callback_without_top_level_stat():
    def callback():
        return 1
    edge=SimpleNamespace(code=callback.__code__,callcount=3,reccallcount=0,inlinetime=.1,totaltime=.2)
    caller=SimpleNamespace(code=test_raw_call_export_retains_callback_without_top_level_stat.__code__,
                           callcount=1,reccallcount=0,inlinetime=.3,totaltime=.5,calls=[edge])
    rows=profiling.raw_profile_rows(SimpleNamespace(getstats=lambda:[caller]))
    assert rows[0]['callees'][0]['code_id']==rows[1]['code_id']
    assert rows[0]['callees'][0]['calls']==3
    assert rows[1]['entry_kind']=='callee_only_no_top_level_stat'
    assert rows[1]['calls'] is None
