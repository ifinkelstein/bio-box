#!/bin/sh
# W-27: fetch the uv.lock wheels of matplotlib and its new dependencies (CPython 3.14, manylinux x86_64),
# check each sha256 against uv.lock, unpack into /tmp/if2453/pydeps_w27. numpy is not fetched (pydeps has it).
set -eu
W=/tmp/if2453/w27/wheels
D=/tmp/if2453/pydeps_w27
mkdir -p "$W" "$D"
while read -r name version url hash; do
  fn=${url##*/}
  [ -f "$W/$fn" ] || curl -fsSL -o "$W/$fn.part" "$url"
  [ -f "$W/$fn" ] || mv "$W/$fn.part" "$W/$fn"
  got=$(sha256sum "$W/$fn" | cut -d' ' -f1)
  want=${hash#sha256:}
  if [ "$got" != "$want" ]; then echo "SHA MISMATCH $fn $got $want"; exit 1; fi
  echo "ok $name $version $got"
  python3 -m zipfile -e "$W/$fn" "$D" 2>/dev/null || unzip -q -o "$W/$fn" -d "$D"
done <<'EOF'
matplotlib 3.11.2 https://files.pythonhosted.org/packages/29/5e/52f56f93d5b20815ee0ae352f7a38eff25580edb180b4fa7df87cfeba68f/matplotlib-3.11.2-cp314-cp314-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl sha256:8c8255de28f986d935a64c9ca71c0ec2d2f41d355691f5ea684725dc91413f71
contourpy 1.4.0 https://files.pythonhosted.org/packages/cf/8c/9677dad226e5aaac5b823d5b09d881ebc32ff1648af417597649841e516c/contourpy-1.4.0-cp314-cp314-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl sha256:e439ab450c93455feb0218532ded3e946ec7a9b5f459069f122cfa44b89229f6
cycler 0.12.1 https://files.pythonhosted.org/packages/e7/05/c19819d5e3d95294a6f5947fb9b9629efb316b96de511b418c53d245aae6/cycler-0.12.1-py3-none-any.whl sha256:85cef7cff222d8644161529808465972e51340599459b8ac3ccbac5a854e0d30
fonttools 4.68.0 https://files.pythonhosted.org/packages/b7/67/e9f2eb2aed41d4ed251a3c6f66e104685f0a03d8ded7ed87f1f0fc63898c/fonttools-4.68.0-cp314-cp314-manylinux2014_x86_64.manylinux_2_17_x86_64.whl sha256:cbb63b26d33acc3eb61f544c22a0cadbe01682fc064b266bb0e33fa66cea853f
kiwisolver 1.5.1 https://files.pythonhosted.org/packages/a1/13/08188f0cafa3a800403e4ff62b9aad4e7a17f9c4c7e080dc8f18c64794cf/kiwisolver-1.5.1-cp314-cp314-manylinux2014_x86_64.manylinux_2_17_x86_64.whl sha256:18a0cfb124546a4c2e6087c5f3029c7f44b37c85b142e0ced71f73a7599ac208
pillow 13.0.0 https://files.pythonhosted.org/packages/1d/f2/1371707ce438ba9a43f9d24089828494ec6f9396e28cd7903aa704dd6306/pillow-13.0.0-cp314-cp314-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl sha256:e91c2b1cbeae768c98752e67e717e283d68dae45717b41f2d91eaee648e37b7d
pyparsing 3.3.3 https://files.pythonhosted.org/packages/38/bb/d215ee7c73b61497b28a5503f9f53523f294fcc936762b7caf90e0c1c2b5/pyparsing-3.3.3-py3-none-any.whl sha256:ece8c00a69cf01b45d0b1dedabb469c90d8caf996d4fda40f147627a122849a4
python-dateutil 4.0.0 https://files.pythonhosted.org/packages/2f/a9/52b89563d154433d79b51ed1eb4842bf7334cb2696891a7fe135618853d7/python_dateutil-4.0.0-py3-none-any.whl sha256:a5e246f79eb80475e691bf93862cead98decd8a498cf53ba0516fcdc0f50dc1c
tzdata 2026.6 https://files.pythonhosted.org/packages/0a/66/9a06c2227fedfd4919f0eaaa8b94377d9730882a7e4f95823c9a331ab5e1/tzdata-2026.6-py2.py3-none-any.whl sha256:0a10eec6e13ef16b627fabed9a5f2ba61a62b1c67f327305f4e7e412a70b4aec
EOF
ls "$D"
