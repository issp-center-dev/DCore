#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Pure-Python checks for the HPhi sector-batching routing used by
``calc_one_body_green_core_parallel`` in internal-loop mode.

These exercise the new grouping / reassembly logic without needing an HPhi
executable: that each task lands in exactly one batch, every same-spin batch is
sector-consistent (one (ex_state, sigma)), cross-spin off-diagonal operators are
singletons, and the task-order reassembly maps every result back correctly.
"""

import itertools

import numpy as np

from dcore.impurity_solvers.hphi_spectrum import _group_tasks_by_sector, _task_key


def _gen_tasks(n_site):
    """Reproduce the task list of calc_one_body_green_core_parallel.gen_p()."""
    n_sigma = 2
    n_excitation = 2
    # p_common contains T_list (a list) -> the full task tuple is unhashable, which
    # is exactly why _task_key drops it; keep it here to mirror the real tasks.
    p_common = (n_site, [1.0], 4, 1e-4, "HPhi", "zvo", "./output", 4)

    def comp(site, sigma):
        return site * n_sigma + sigma

    tasks = []
    for si, sgi in itertools.product(range(n_site), range(n_sigma)):
        for sj, sgj in itertools.product(range(n_site), range(n_sigma)):
            if comp(si, sgi) > comp(sj, sgj):
                continue
            for idx_flg in range(2):
                for ex_state in range(n_excitation):
                    if comp(si, sgi) == comp(sj, sgj) and idx_flg == 1:
                        continue
                    tasks.append((si, sgi, sj, sgj, idx_flg, ex_state, p_common))
    return tasks


def test_group_tasks_cover_each_task_exactly_once():
    tasks = _gen_tasks(2)
    batches = _group_tasks_by_sector(tasks)
    flat = [t for b in batches for t in b]
    assert len(flat) == len(tasks)
    assert {_task_key(t) for t in flat} == {_task_key(t) for t in tasks}
    # no task duplicated across or within batches
    assert len({_task_key(t) for t in flat}) == len(flat)


def test_batches_are_sector_consistent():
    tasks = _gen_tasks(2)
    for b in _group_tasks_by_sector(_gen_tasks(2)):
        sector_keys = set()
        for (si, sgi, sj, sgj, idx_flg, ex_state, _) in b:
            if sgi == sgj:
                sector_keys.add((ex_state, sgi))
            else:
                # cross-spin operators mix Sz -> must be alone in their launch
                assert len(b) == 1, "cross-spin operator was batched with others"
        # a same-spin batch carries exactly one (ex_state, sigma)
        assert len(sector_keys) <= 1


def test_task_order_reassembly_round_trips():
    tasks = _gen_tasks(2)
    batches = _group_tasks_by_sector(tasks)

    def synth(t):
        return np.array(_task_key(t), dtype=float)

    # emulate calc_one_body_green_batch: each batch returns {task_key: result}
    result_map = {}
    for b in batches:
        for t in b:
            result_map[_task_key(t)] = synth(t)

    assert len(result_map) == len(tasks)
    results = [result_map[_task_key(t)] for t in tasks]
    for t, r in zip(tasks, results):
        assert np.array_equal(r, synth(t))
