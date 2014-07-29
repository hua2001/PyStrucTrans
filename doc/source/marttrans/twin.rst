Twin
====

Martensite variants form twins.
The module :py:mod:`pystructrans.marttrans.twin` contains two main classes

* class :py:class:`pystructrans.TwinSystem`
* class :py:class:`pystructrans.marttrans.twin.TwinPair`

TwinPair object
---------------
The class :py:class:`pystructrans.marttrans.twin.TwinPair`

.. autoclass:: pystructrans.marttrans.twin.TwinPair
    :members: getUi, getUj, iscompound, istypeI, istypeII, iscompatible

    :param Ui: one martensite variant
    :param Uj: another martensite variant
    :param boolean skipcheck: skip the input checking, default is ``False``. Use with caution!
    :raises ValueError: if Ui or Uj is not a 3 x 3 positive definite symmetric matrix.

    .. automethod:: pystructrans.marttrans.twin.TwinPair.istwinnable

        The variants :math:`\mathbf U_i` and :math:`\mathbf U_j` are twinnable if

        #. :math:`\mathbf U_i \neq \mathbf U_j`
        #. The matrix :math:`\mathbf C = \mathbf U_i^{-1}\mathbf U_j^2\mathbf U_i^{-1}`
           has the middle eigenvalue **1**.

        :rtype: boolean

    .. automethod:: pystructrans.marttrans.twin.TwinPair.gettwinparam

        Twin parameters are the :math:`\mathbf a` and :math:`\mathbf n`
        in the two solutions to the the twinning equation:

        .. math::

            \mathbf Q\mathbf U_i - \mathbf U_j = \mathbf a \otimes \mathbf n

        For a compound twin, these two sets of solutions degenerate

        :return: ((:math:`\mathbf Q_\text{I}`, :math:`\mathbf a_\text{I}`, :math:`\mathbf n_\text{I}`),
                 (:math:`\mathbf Q_\text{II}`, :math:`\mathbf a_\text{II}`, :math:`\mathbf n_\text{II}`))

        :raises AttributeError: if not twinnable

    .. automethod:: pystructrans.marttrans.twin.TwinPair.habitplanes

        :param str twintype: "I", "II" or "C". If provided only associated volume fractions will be returned.
        :return: If ``twintype`` is "C", return
                 (((:math:`\mathbf R_\text{I1}`, :math:`\mathbf b_\text{I1}`, :math:`\mathbf m_\text{I1}`),
                 (:math:`\mathbf R_\text{I2}`, :math:`\mathbf b_\text{I2}`, :math:`\mathbf m_\text{I2}`)),
                 ((:math:`\mathbf R_\text{II1}`, :math:`\mathbf b_\text{II1}`, :math:`\mathbf m_\text{II1}`)
                 (:math:`\mathbf R_\text{II2}`, :math:`\mathbf b_\text{II2}`, :math:`\mathbf m_\text{II2}`)))
                 otherwise return only the first two or the last two.
        :raises AttributeError: if not compatible
        :raises ValueError: unknown ``twintype``

    .. automethod:: pystructrans.marttrans.twin.TwinPair.isconventional

        A twin pair is conventional if the 180 |deg| rotation relating the two variant
        is in the Laue group of the high symmetry (parent) phase.

        :param lauegroup: the Laue group. If not given, use full cubic Laue group
        :rtype: boolean

    .. automethod:: pystructrans.marttrans.twin.TwinPair.volfrac

        :param str twintype: "I", "II" or "C". If provided only associated volume fractions will be returned.
        :return: If ``twintype`` is "C", return ((`f` :sub:`I`, 1 - `f` :sub:`I`), (`f` :sub:`II`, 1 - `f` :sub:`II`)),
                 otherwise return only (`f` :sub:`I`, 1 - `f` :sub:`I`) or (`f` :sub:`II`, 1 - `f` :sub:`II`).

    .. automethod:: pystructrans.marttrans.twin.TwinPair.XI

        Let :math:`\mathbf t` be the axis of 180 |deg| rotation relating the two variants.
        The cofactor parameter :math:`X_\text{I} = \left|\mathbf U_i^{-1}\mathbf t\right| = \left|\mathbf U_j^{-1}\mathbf t\right|`.

        :raise AttributeError: if is compound twin

    .. automethod:: pystructrans.marttrans.twin.TwinPair.XII

        Let :math:`\mathbf t` be the axis of 180 |deg| rotation relating the two variants.
        The cofactor parameter :math:`X_\text{II} = \left|\mathbf U_i\mathbf t\right| = \left|\mathbf U_j\mathbf t\right|`.

        :raises AttributeError: if is compound twin

.. |deg| unicode:: 0xB0 .. degree sign