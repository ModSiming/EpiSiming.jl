# EpiSiming.jl

[docs-dev-img]: https://img.shields.io/badge/docs-dev-green.svg
[docs-dev-url]: https://modsiming.github.io/EpiSiming.jl/dev/


[![][docs-dev-img]][docs-dev-url] [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-orange.svg)](https://www.gnu.org/licenses/gpl-3.0) ![GitHub repo size](https://img.shields.io/github/repo-size/modsiming/episiming.jl)

Epidemics simulator (Julia version).

This is the Julia version of [ModSiming/EpiSiming](https://github.com/ModSiming/EpiSiming).

## Team

- [Beatriz Farah](https://github.com/beafarah)
- [Cynthia Herkenhoff](https://github.com/herkenhoff-cynthia)
- [Gil Miranda](https://github.com/mirandagil/)
- [Ricardo Rosa](http://www.im.ufrj.br/rrosa/)
- [Rodrigo Peregrino](https://github.com/rodlcp)
- [Thiago Holleben](https://github.com/hollebenthiago)

## About the name EpiSiming

The name **EpiSiming** comes from combining the prefix **Epi**, from epidemic (of course), with the chinese word **Siming**, which, [according to Wikipedia](https://en.wikipedia.org/wiki/Siming_(deity)), refers to a Chinese deity that makes fine adjustments to the fate of humankind.

*Siming* is also the name of a Chinise district, and the pronunciation of this Chinese district can be heard in [How to Pronounce Siming District - PronounceNames.com](https://www.youtube.com/watch?v=VXkclmg96BQ).

## A random scenario example

The core engine is implemented, but so far we have only build one scenario, which is a random scenario on a rectangular area. A scenario for the city of Rio de Janeiro has been implemented in the python version of the project and will soon be adapted to julia and included here.

For the random scenario, have a look [Epidemic simulation via a discrete-time, agent-based stochastic model with a random scenario](docs/src/examples/random_scenario.md).

## License

EpiSiming.jl - Epidemics simulator

Copyright (C) 2021 by Ricardo Martins da Silva Rosa (IM/UFRJ) and other members of the github organization ModSiming/EpiSiming.jl

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
