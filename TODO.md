# ToDo

- [ ] Port all heavy numerical operations (especially in `Spicy.Math`) to `runQ` enabled Accelerate
- [ ] Define types for QC wrappers, that are easily extensible by configuration files and can be parsed by YAML files at run time
- [ ] Rewrite the `Spicy.MolecularSystem` module and split it into submodules. Use better naming scheme for functions.
- [ ] Write wrappers for quantum chemistry code
- [ ] Define a YAML based input format
- [ ] Implement ONIOM methods
  - [ ] ONIOM2
    - [ ] Gradients
  - [ ] ONIOM3
    - [ ] Gradients
  - [ ] Arbitrary high ONIOM
  - [ ] Multicentre ONIOM
  - [ ] Excited state ONIOM
  - [ ] Dynamic MD ONIOM
- [ ] Implement Quantum Dynamics MD algorithm
- [ ] (Re)Write tests
- [ ] (Re)Write benchmarks
- [ ] Rewrite parsers and writers
