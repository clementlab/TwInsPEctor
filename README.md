# TwInsPEctor

Tools and scripts for analyzing twin prime editing experiments, including data processing and visualization.

## Features

- Command-line interface for streamlined analysis
- Data processing utilities for twin prime editing experiments
- Visualization tools for editing outcomes

## Installation

### Using Conda

```bash
conda install bioconda::twinspector
```

### From Source

```bash
git clone https://github.com/clementlab/TwInsPEctor.git
cd TwInsPEctor
pip install .
```

## Usage

After installation, use the CLI:

```bash
TwInsPEctor --help
```

Or run the main module directly:

```bash
python -m TwInsPEctor
```

## Requirements

- Python >=3.8
- [CRISPResso2](https://github.com/pinellolab/CRISPResso2) (installed via Bioconda)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Authors

- Nate Masson
- Kendell Clement

## Dependency Notice

This software requires [CRISPResso2](https://github.com/pinellolab/CRISPResso2) to be installed separately.

CRISPResso2 is distributed under its own license terms, which may
restrict commercial use.

Users are responsible for ensuring compliance with the CRISPResso2
license when using this software.

This project does not redistribute CRISPResso2 and does not grant
any rights to it.
