# KDUQ-TWOFNR
Attemting to run KDUQ parameters into TWOFNR for ADWA, then through either TWOFMR or maybe DWUCK4 for use in 19O(d,p) data analysis. 


### Outlining:
```mermaid
graph TB

Generate_KDUQ_Samples --> Build_FRONT_program --> Run_FRONT_program -- no V&F --> Run_TWOFNR_program
Run_FRONT_program -- with V&F --> Convert_to_DWUCK4 --> Run_DWUCK4
Run_TWOFNR_program --> Plot_all_outputs
Run_DWUCK4 --> Plot_all_outputs

```

working in snakefile:
- [x] Generate KDUQ samples
- [ ] build front program
- [ ] run front program
- [ ] run twofnr program
	- [ ] OR!
		- [ ] convert to dwuck
		- [ ] run dwuck4
- [ ] plot all outputs


