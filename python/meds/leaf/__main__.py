"""`python -m meds.leaf` — run the round-trip self-test against libmeds_plant_c."""
from . import self_test

flux = self_test()
print(f"meds.leaf self-test PASS: A_net={flux.A_net:.3f} gs={flux.gs:.4f} "
      f"Ci={flux.ci:.2f} limitation={flux.limitation.name} converged={flux.converged}")
