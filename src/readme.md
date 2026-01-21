## Overview of current and future changes
This document focuses on the C++ rewrite for input, mapping (seed, chain, align) and output. Building the index is not included yet.

### Background
Most of the C code to change is found in `map.cpp` and its dependencies. Importantly, we are concerned with `mm_map_file_frag`, `kt_pipeline`/`worker_pipeline`, and `kt_for`/`worker_for`, and the functions they call.
- `mm_map_file_frag`: uses `kt_pipeline` to call `worker_pipeline` after some initial setup
- `worker_pipeline`: 3 main steps, input, mapping and output. Uses `kt_for` to call `worker_for`.
- `worker_for`: heavy computation steps: where seed, chain and align happens per segment/fragment

Currently, mapping is optimized for CPU multithreading where we have a thread pool that completes each segment/fragment one at a time, reusing thread-local buffers. For example, thread1 (t1) does seed, chain, align on segment1, t2 on segment2, t3 on segment3. Then, t1 finishes and moves on to segment4. 

Instead, we want threads to finish step by step, i.e. finish all of seeding before moving to chaining. This is important for GPU batching since data parallelism exists within each of seed, chain and align steps. However, it is a problem in multithreaded implementations since each step generates some intermediate data that is backed by each thread's own memory pool. Hence, we want to introduce the idea of batches where t1 owns segment [0, n), t2 owns [n, 2n) and so on. 

### Design overview
We separate the algorithms and the data they operate on. Algorithm-heavy classes include the seeder, chainer, aligner and mapper, and the data is put into custom types passed between those classes. We also want to rewrite the input and output to ensure data is in memory as needed.

#### `FileReader`
`FileReader` is the main input reader for query files and it produces the expected input for algorithm classes
```c++
// list of input files (FASTA/FASTQ) to read 
std::vector<std::string> files = {"f1.fa", "f2.fa"};
// set config
FileReaderConfig config = {
    .enable_quality = opt->flag & MM_F_OUT_SAM && !(opt->flag & MM_F_NO_QUAL),
    .enable_comment = opt->flag & MM_F_COPY_COMMENT,
    .fragment_mode = opt->flag & MM_F_FRAG_MODE,
};
std::unique_ptr<FileReader> file_reader(files, config);
// read the next n segments, allowing to break up the input into multiple segments for multithreading
shared_ptr<InputDataFragments> input = file_reader->readNextSegments(n);
```
A single instance is intended to be created before `worker_pipeline` runs, and each batch will read some number of segments in step 0 of `worker_pipeline` to create a `MappingContext` (see below) for next steps

#### `Mapper`
`Mapper` is the main interface intended to replace `worker_for`, calling the seeder, chainer and aligner. It operates on `MappingContext` which contains the input reads, configs, and the index. An example usage is:
```c++
// create configs
MapperConfig cfg = {
    .flags = p->opt->flag,
    .paired_end_orientation = bitset<2>(opt->pe_ori),
    .seed_cfg = MapperConfig::SeederConfig{
        .query_occurrence_fraction = opt->q_occ_frac,
        .seed_occurrence_threshold = opt->mid_occ,
        .hard_seed_occurrence_threshold = opt->max_max_occ,
        .seed_occurrence_distance = opt->occ_dist,
        .sdust_threshold = opt->sdust_thres,
    },
    .chain_cfg = MapperConfig::ChainerConfig{
        .max_query_gap = opt->max_gap,
        .max_ref_gap = opt->max_gap_ref,
        .max_fragment_length = opt->max_frag_len,
        .max_skip = opt->max_chain_skip,
        .max_predecessors = opt->max_chain_iter,
        .bandwidth = opt->bw,
        .bandwidth_long = opt->bw_long,
        .min_chain_anchors = opt->min_cnt,
        .min_chain_score = opt->min_chain_score,
        .chain_gap_scale = opt->chain_gap_scale,
        .chain_skip_scale = opt->chain_skip_scale,
    }
    // TODO: .align_cfg = {},
};
// create a context object: obj that has config and input to run seed,chain,align on
auto context = make_shared<MappingContext>(
    cfg,
    shared_ptr<mm_idx_t>(const_cast<mm_idx_t *>(p->mi), [](mm_idx_t *) {}), // index built in index.c
    input); // input object from FileReader
// Create a mapper object
Mapper mapper();
mapper->map(context); // this does seed, chain and align for this context object
```

#### `Seeder`, `Chainer`, `Aligner`
These classes are meant to be called by `Mapper::map`. The idea is that the context object passed in will be completed in batches, hence there needs to be a batching step prior to the `Seeder` running but is not yet done. See TODOs in `Mapper::map` for more info on the flow of data.

#### Misc. files
We also have `types.hpp`, `types.cpp` and `utils.hpp` that provide some support for all the above classes

### Current progress
**Legend**
- ⬜ Not started
- 🟨 In progress / partially complete
- ✅ Done but not tested
- 🧪 Done and tested
------------------
- **FileReader** 🟨  
    - FASTA support: 🧪  
    - FASTQ support: ✅
    - Compressed files (e.g. `.gz`): ⬜  
- **Seeder** 🧪  
- **Chainer** 🟨  
    - DP-based score population: 🧪  
    - Backtracking & chain generation: 🧪  
    - RMQ-based score population: ⬜
    - Redo chaining: ⬜
- **Aligner** ⬜  
- **Mapper** 🟨  
    - Batching ⬜  
    - Algorithm classes: 🟨
    - Helper functions: 🧪  
    - Return type: ⬜
- **Output** ⬜

 There still exists functions that haven't been properly designed/thought of. For example, `merge_hits` which is in step 1 of `worker_pipeline`. 
 
 ### Future works
 - Merge this batched CPU processing with GPU chaining and seeding
 - Rewrite index.c