// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

// Some tests (scripts, the README example) change the working directory of the process, which is
// shared by every thread: the test classes run one after the other.
[assembly: CollectionBehavior(DisableTestParallelization = true)]
