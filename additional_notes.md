# Additional Notes

## Input file patterns

In `config.toml`, the comment for `fpath` is incorrect:

```toml
fpath = "/path/to/model/output"   # must contain cobalt_btm.nc, cobalt_tr.nc, mom6.nc
```

The correct input file patterns are `*cobalt_btm.nc`, `*cobalt_tracers.nc`, and `*ocean_daily.nc`, as hard-coded in `model_reader.py` ([lines 199–215](https://github.com/dksasaki/cbed_R2py/blob/main/scripts/model_reader.py#L199-L215)). The variables read from each file are also hard-coded ([lines 99–143](https://github.com/dksasaki/cbed_R2py/blob/main/scripts/model_reader.py#L99-L143)).

Each file must have the following dimensions:

- `*cobalt_btm.nc` : `(time, yh, xh)`
- `*cobalt_tracers.nc` : `(time, z_l, yh, xh)` — bottom-most level is extracted
- `*ocean_daily.nc` : `(time, zl, yh, xh)` — bottom-most level is extracted



## Time averaging inconsistency

There is an inconsistency in where the time average is computed. In `model_reader.py`, the average is only applied inside `_read_mom`, while the averages for the COBALT datasets are deferred to `cbed_wrapper.py` ([lines 258–261](https://github.com/dksasaki/cbed_R2py/blob/main/scripts/cbed_wrapper.py#L258-L261)):

```python
dscob  = ds_dict['dscobalt_btm'].mean(dim='time')
dscob2 = ds_dict['dscobalt_tr'].mean(dim='time')
dsmom  = ds_dict['dsmom']
```

Ideally, all time averages should be computed inside `model_reader.py` at read time, alongside `_read_mom`.

## Vertical fill

Considering `*cobalt_tracers.nc` , before selecting the bottom-most level we [extrapolate the last valid result in cobalt fields to the bottom most layer](https://github.com/dksasaki/cbed_R2py/blob/main/scripts/model_reader.py#L213-L214). This propagates valid values into NaN-padded layers, which is necessary when the model's land mask leaves trailing NaNs at depth — without it, `isel(z_l=-1)` could return NaN for valid ocean points.


## other comments
### in [README](https://github.com/dksasaki/cbed_R2py/blob/main/README.md)

Section 5. Scripts
- `chunk=0` runs the full domain

	you shoud use a batch job to run separate instances of the following command
	```shell
	python cbed_wrapper.py 200 2  # uses 200 processors and selects the second chunk
	python cbed_wrapper.py 200 3  # uses 200 processors and selects the third chunk
	```

### Global variables in worker processes

The datasets `dsmom_a`, `dscob_a`, `dscob2_a`, `dsporo_a`, and `valid_points` are [accessed as implicit globals](https://github.com/dksasaki/cbed_R2py/blob/main/scripts/cbed_wrapper.py#L192-#L200) inside  `run_point`. This works because `multiprocessing` forks the parent process, **but it is fragile** — any refactoring that moves `run_point` out of the main script would require making these dependencies explicit.

### Comment on the wrapper
We are using *literally* r within python. [Each r task runs about 20 tasks and then resets to avoid overhead issues.](https://github.com/dksasaki/cbed_R2py/blob/main/scripts/cbed_wrapper.py#L307-L309) This is  needed to control memory growth from repeated R/rpy2 calls within a single process. 
- `maxtasksperchild=20`

	
