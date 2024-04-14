# RM-17O-by-VCOF-CRDS

Data and code repository for *Linking the oxygen-17 compositions of water and carbonate reference materials using infrared absorption spectroscopy of carbon dioxide* (Chaillot et al., in review [10.31223/X53Q43](https://doi.org/10.31223/X53Q43)).

## Installation and execution

* Install the [pixi](https://pixi.sh) package management tool.
* To run the whole code from scratch: `pixi run all`
* To run specific parts of the code, use the following commands.

```sh
pixi run stability-wg
pixi run memory-effect
pixi run pressure-effect
pixi run linearity
pixi run d13C-effect
pixi run carbonate-repeatability
pixi run ref-materials
pixi run pub-comparison  # reads the output of ref-materials
pixi run all-together    # reads the output of linearity and pub-comparison
pixi run stdz-example    # reads the output of linearity
```

## Contents


| Direcory                  | Content                                                      |
|:------------------------- |:-------------------------------------------------------------|
|`0-stdz_example`           | Toy example illustrating our standardization methods         |
|`1-stability-wg`           | Characterize instrumental stability                          |
|`2-memory-effect`          | Test for memory effects                                      |
|`3-pressure-effect`        | Test for pressure effects                                    |
|`4-linearity`              | Test Δ17O linearity based on water equilibration experiments |
|`5-d13C-effect`            | Check whether δ13C affects Δ17O measurements                 |
|`6-carbonate-repeatability`| Quantify acid bath repeatability                             |
|`7-ref-materials`          | Characterize water & carbonate ref materials                 |
|`8-pub-comparison`         | Compare our results with previous publications               |
|`9-all-together`           | Compile additional plots and tables                          |
|`lib`                      | Shared libraries                                             |


## Contact

Please feel free to open an issue on [GitHub](https://github.com/mdaeron/RM-17O-by-VCOF-CRDS) or contact [M. Daëron](daeron@lsce.ipsl.fr) directly.
