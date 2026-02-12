# Results and validation

This document presents validation results for fiber generation using two methods: the Bayer method and the Doste method.

The `examples` folder contains two geometries: a truncated biventricle (`examples/truncated`) for the Bayer method and a biventricle with outflow tracts (`examples/ot`) for the Doste method. 

Validation scripts (`validation_bayer.py` and `validation_doste.py`) can be run to generate the correlation plots shown below. ParaView Python scripts (`paraview_bayer.py` and `paraview_doste.py`) can be used to generate the different visualization views (fiber, sheet, sheet-normal orientations) for each case.

The validation codes first calculate the $\alpha$ and $\beta$ angles using scalar interpolations. Then, the $\alpha$ and $\beta$ angles are calculated from the fiber direction $\mathbf f$ and the orthogonal basis $(\mathbf e_c, \mathbf e_\ell, \mathbf e_t)$. 

We do this for three combinations of parameters:
1. Setting all $\alpha$ to default values, and all $\beta$ to 0.
2. Setting all $\alpha$ to 0, and all $\beta$ to default values.
3. Setting all $\alpha$ and $\beta$ to default values.

This allows to isolate and identify issues in the $\alpha$ and $\beta$ rotations. 

Ideally both scalar interpolated and fiber derived angles should match exactly, but given that the orthogonal basis $(\mathbf e_c, \mathbf e_\ell, \mathbf e_t)$ also needs to be interpolated some differences arise. 


Notes:
- To run the `validation*.py` scripts the file with the Laplace results must be created. This can be done using the `main*.py` scripts.
- To run the `paraview*.py` you must run the `validation*.py` codes first. 
- `paraview*.py` codes are run within the Paraview GUI. To do so, 
   1. Open Paraview, go to the Python Shell (if not visible, go to View, and click so it appears, usually at the bottom panel).
   2. Click on Run Script and select the desire script. 

---

## Bayer Method

The Bayer method results are demonstrated on a truncated biventricular geometry.

### Fiber Orientation

![Bayer Fiber Full View](example/biv_truncated/bayer_fiber.png)

*Figure 1: Fiber orientation field generated using the Bayer method - full view*

![Bayer Fiber Slice View](example/biv_truncated/bayer_fiber_slice.png)

*Figure 2: Fiber orientation field generated using the Bayer method - slice view*

### Sheet Orientation

![Bayer Sheet Full View](example/biv_truncated/bayer_sheet.png)

*Figure 3: Sheet orientation field generated using the Bayer method - full view*

![Bayer Sheet Slice View](example/biv_truncated/bayer_sheet_slice.png)

*Figure 4: Sheet orientation field generated using the Bayer method - slice view*

### Sheet-Normal Orientation

![Bayer Sheet-Normal Full View](example/biv_truncated/bayer_sheet-normal.png)

*Figure 5: Sheet-normal orientation field generated using the Bayer method - full view*

![Bayer Sheet-Normal Slice View](example/biv_truncated/bayer_sheet-normal_slice.png)

*Figure 6: Sheet-normal orientation field generated using the Bayer method - slice view*

### Angle Correlations

To check the code, we first calculate the $\alpha$ and $\beta$ angles using scalar interpolations. Then, we calculate the $\alpha$ and $\beta$ angles using the fiber direction $\mathbf f$ and the orthogonal basis $\mathbf e_c$, $\mathbf e_\ell$, $\mathbf e_t$. 

![Bayer Angle Correlations](example/biv_truncated/bayer_angle_correlations.png)

*Figure 7: Correlation plots comparing scalar interpolation angles with fiber derived angles for the Bayer method. Blue and red dots show the $\alpha$ and $\beta$ angles. For reference, the original Bayer method with no modifications is shown.*

---

## Doste Method

The Doste method results are demonstrated on a complete biventricular geometry with outflow tracts.

### Fiber Orientation

![Doste Fiber Full View](example/biv_with_outflow_tracts/doste_fiber.png)

*Figure 8: Fiber orientation field generated using the Doste method - full view*

![Doste Fiber Slice View](example/biv_with_outflow_tracts/doste_fiber_slice.png)

*Figure 9: Fiber orientation field generated using the Doste method - slice view*

### Sheet Orientation

![Doste Sheet Full View](example/biv_with_outflow_tracts/doste_sheet.png)

*Figure 10: Sheet orientation field generated using the Doste method - full view*

![Doste Sheet Slice View](example/biv_with_outflow_tracts/doste_sheet_slice.png)

*Figure 11: Sheet orientation field generated using the Doste method - slice view*

### Sheet-Normal Orientation

![Doste Sheet-Normal Full View](example/biv_with_outflow_tracts/doste_sheet-normal.png)

*Figure 12: Sheet-normal orientation field generated using the Doste method - full view*

![Doste Sheet-Normal Slice View](example/biv_with_outflow_tracts/doste_sheet-normal_slice.png)

*Figure 13: Sheet-normal orientation field generated using the Doste method - slice view*

### Angle Correlations

![Doste Angle Correlations](example/biv_with_outflow_tracts/doste_angle_correlations.png)

*Figure 14: Correlation plots comparing scalar interpolation angles with fiber derived angles for the Doste method. Blue and red dots show the $\alpha$ and $\beta$ angles**

---