Processes, formulations and numerics
====================================

Transport equation
------------------

.. math::
   q = \alpha C_b \frac{\rho}{g} \sqrt{\frac{d}{D}} \left ( u - u_{th} \right ) ^3

   \alpha = \left ( \frac{0.174}{log \left ( \frac{z_0}{k} \right )} \right ) ^3

Velocity thresholds
-------------------

Grain size
^^^^^^^^^^

Bed slope
^^^^^^^^^

Moisture content
^^^^^^^^^^^^^^^^

Multiple sediment fractions
---------------------------

Sediment transport of a multiple fractions is computed equal to transport of a single fraction.
For each fraction transport is computed given an equilibrium concentration :math:`C_{eq,i}`.
Since not all fractions can be at equilibrium concentration at the same time at the same place, a distribution is needed that favors one fraction over another.
This is the key to multiple sediment fraction transport.

A common approach is to use the sediment distribution of the bed to weigh the equilibrium concentrations of the different fractions [Delft3DManual]_.
This approach has the disadvantage that not only pickup of sediment, but also deposition of sediment depends on the sediment distribution in the bed.
This is physically incorrect since deposition is not strongly correlated to the bed composition.

A concequence of this approach is that fine sediments picked up downwind cannot pass a patch of coarse sediment where the fraction of fine sediments is small.
In this patch, the equilibrium concentration of fines is artificially lowered and fines are deposited, while coarse sediment is picked up.

To overcome this behavior the distribution for erosion and deposition is determined separately according to equation :eq:`distribution`.
First the distribution in the air :math:`p_{air,i}` is determined.
This distribution can sum up to a value smaller than, equal to or larger than unity that correspond to three regimes (:math:`N` is the number of fractions):

* :math:`\sum_{i=0}^N{p_{air,i}} < 1.0`: less than capacity transport, erosion may occur
* :math:`\sum_{i=0}^N{p_{air,i}} = 1.0`: capacity transport, no exchange with the bed
* :math:`\sum_{i=0}^N{p_{air,i}} > 1.0`: more than capacity transport, deposition will occur

In case of deposition the distribution in the air is normalized and used as distribution in the transport formulation.
In case of erosion the distribution in the bed :math:`p_{bed,i}` is used to *fill* the available space in the air column.
The distribution of the bed is determined by the mass of each sediment fraction :math:`m_i` available in the top layer of the bed.
This means that if a grid cell is at 70% of capacity transport, 30% of the distribution used in the transport formulation is determined by the distribition in the bed and 70% is determined by the distribution in the air.

.. math::
   :label: distribution

   p_i = \left \{  \begin{array}{l l}
      \frac{p_{air,i}}{\sum_{i=0}^N{p_{air,i}}} &
         \quad \text{if } \sum_{i=0}^N{p_{air,i}} \geq 1.0 \text{ (deposition)} \\
      p_{air,i} + p_{bed,i} \left (1 - \sum_{i=0}^N{p_{air,i}} \right ) &
         \quad \text{if } \sum_{i=0}^N{p_{air,i}} < 1.0 \text{ (erosion)}
      \end{array} \right.

   p_{air,i} = \frac{C_{t,i}}{C_{u,i}} \quad ; \quad
   p_{bed,i} = \frac{m_i}{\sum_{i=0}^N{m_i}}

   
Hydraulic mixing and depostion
------------------------------

The left-hand boundary can be taken as a sea boundary imposing tides and waves onto the sandy bed.
Apart from sheltering the bed from the wind and thus transport, it may also simulate processes of mixing and deposition.

If grid cells are flooded the top layers of the bed are mixed to simulate stirring by waves.
Mixing takes place over a certain depth indicated by the depth of disturbance (DOD).
Over the depth of disturbance the average sediment distribution is determined and the actual sediment distributions are then replaced by its average.
The depth of disturbance is computed based on the offshore wave height :math:`H_s`, but it is maximized by a maximum wave height over depth ratio :math:`\gamma` according to equation :eq:`DOD`. The factor :math:`f_{DOD}` is an empirical ratio between onshore wave height and depth of disturbance which is 10% - 20% [Masselink2007]_. :math:`\eta` and :math:`z` are the water bed level respectively.

.. math::
   DOD(t) = f_{DOD} \cdot \min(H_s, \gamma (\eta(t) - z(t)))
   :label: DOD

Flooding of grid cells may also lead to a fresh sediment deposit at these grid cells.
The deposit is computed as adaptation of the supply term :math:`S` according to equation :eq:`deposit`.
The deposit is computed by multiplying the fall velocity :math:`w` with the sediment concentration in the water column :math:`C_w` and the time step :math:`\Delta t`.
Deposition is maximized by the water depth :math:`\eta - z` and weighted according to the initial sediment distribution of the bed :math:`p_{i,initial}` over :math:`N` fractions.
           
.. math::
   S_i(t) = S_i(t) - \frac{C_w \cdot p_{i,initial}}{\sum_{i=0}^N{p_{i,initial}}} \min( w \Delta t ; \eta(t) - z(t) )
   :label: deposit

References
----------

.. [Delft3DManual] Delft3D manual
.. [Masselink2007] Paper about DOD
