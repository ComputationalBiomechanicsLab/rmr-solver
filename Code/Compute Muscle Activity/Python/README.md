# Installation

- Install python3.8 (required for OpenSim 4.3, see: https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+in+Python)
- Open terminal in this folder
- (*optional*) Create a virtual environment `python -m venv venv/`. Activate with `source venv/Scripts/activate`
- Install pip packages with `pip install -r requirements_new.txt`
- Install OpenSim 4.3 python module with:

```bash
# install python module with pip
pip install /c/OpenSim\ 4.3/sdk/Python
```
- Test with `python -c 'import opensim' && echo "It works!"
