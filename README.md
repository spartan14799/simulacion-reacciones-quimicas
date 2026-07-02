
## run

install uv and pixi
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
powershell -c "irm https://pixi.sh/install.ps1 | iex"
```

instalar paquete editable dentro de pixi
```powershell
pixi run pip install -e .
```

close and run in new terminal
```powershell
pixi run jupyter notebook notebooks/<notebook>
```