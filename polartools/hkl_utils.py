"""
Auxilary HKL functions.

.. autosummary::
    ~sampleChange
    ~sampleList
    ~list_reflections
    ~or_swap
    ~setor0
    ~setor1
    ~set_orienting
    ~list_orienting
    ~or0
    ~or1
    ~compute_UB
    ~calc_UB
    ~setmode
    ~ca
    ~br
    ~uan
    ~wh
    ~setlat
"""

import bluesky.plan_stubs as bps

try:
    import gi

    gi.require_version("Hkl", "5.0")
    from hkl import cahkl
    from hkl.user import _check_geom_selected, _geom_
except ModuleNotFoundError:
    print("gi module is not installed, the hkl_utils functions will not work!")
    cahkl = _check_geom_selected = _geom_ = None


"""
Most of the functions below are only working for the six circle diffractometer
(diffract) right now. This will be changed ...
"""


def sampleChange(sample_key=None):
    """
    Change selected sample in hklpy.

    Parameters
    ----------
    sample_key : string, optional
        Name of the sample as set in hklpy. If None it will ask for which
        sample.
    """

    if sample_key is None:
        d = _geom_.calc._samples.keys()
        print("Sample keys:", list(d))
        sample_key = input(
            "\nEnter sample key [{}]: ".format(_geom_.calc.sample.name)
        )
        if not sample_key:
            sample_key = _geom_.calc.sample.name
    try:
        _geom_.calc.sample = _geom_.calc._samples[
            sample_key
        ]  # define the current sample
        print("\nCurrent sample: " + _geom_.calc.sample.name)
        compute_UB()

    except KeyError:
        print("Not a valid sample key")


def sampleList():
    """List all samples currently defined in hklpy; specify  current one."""
    samples = _geom_.calc._samples
    print("")
    for x in list(samples.keys())[1:]:
        orienting_refl = samples[x]._orientation_reflections
        print("Sample = {}".format(x))
        print("Lattice: ", end="")
        [print(num, end=" ") for num in samples[x].lattice._fields]
        print("=", end=" ")
        [print(num, end=" ") for num in samples[x].lattice]
        print("")
        for i, ref in enumerate(samples[x]._sample.reflections_get()):
            if orienting_refl[0] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>3}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   ".format(
                        "or0",
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                    )
                )
            elif orienting_refl[1] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>3}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}  ".format(
                        "or1",
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                    )
                )
        print(
            "======================================================================"
        )
    print("\nCurrent sample: " + _geom_.calc.sample.name)


def list_reflections(all_samples=False):
    """
    Lists all reflections in defined in hklpy.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    all_samples : boolean, optional
        If True, it will list the reflections for all samples, if False, only
        the current sample. Defaults to False.
    """
    _check_geom_selected()
    if all_samples:
        samples = _geom_.calc._samples.values()
    else:
        samples = [_geom_.calc._sample]
    for sample in samples:
        print("Sample: {}".format(sample.name))
        orienting_refl = sample._orientation_reflections
        print(
            "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
                "#",
                "H",
                "K",
                "L",
                "Delta",
                "Theta",
                "Chi",
                "Phi",
                "Gamma",
                "Mu",
                "orienting",
            )
        )
        for i, ref in enumerate(sample._sample.reflections_get()):
            if orienting_refl[0] == ref:
                # or0_old = i
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                        "first",
                    )
                )
            elif orienting_refl[1] == ref:
                # or1_old = i
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                        "second",
                    )
                )
            else:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                    )
                )
        if len(samples) > 1 and all_samples:
            print(
                "============================================================================"
            )


def or_swap():
    """Swaps the two orientation reflections in hklpy."""
    sample = _geom_.calc._sample
    sample.swap_orientation_reflections()
    list_reflections()
    print("Computing UB!")
    sample.compute_UB(
        sample._orientation_reflections[0], sample._orientation_reflections[1]
    )
    _geom_.forward(1, 0, 0)


def setor0(*args):

    """
    Sets the primary orientation in hklpy.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    delta, th, chi, phi, gamma, mu : float, optional
        Values of motor positions for current reflection. If None, it will ask
        for it.
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """

    sample = _geom_.calc._sample
    ref = sample._sample.reflections_get()
    orienting_refl = sample._orientation_reflections

    if _geom_.name == "diffract" and len(args) == 9:
        h, k, l, delta, th, chi, phi, gamma, mu = args
    else:
        if len(orienting_refl) > 1:
            for ref in sample._sample.reflections_get():
                if ref == orienting_refl[0]:
                    pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                    old_delta = pos[4]
                    old_th = pos[1]
                    old_chi = pos[2]
                    old_phi = pos[3]
                    old_gamma = pos[5]
                    old_mu = pos[0]
                    old_h, old_k, old_l = ref.hkl_get()

        else:
            old_delta = 60
            old_th = 30
            old_chi = 0
            old_phi = 0
            old_gamma = 0
            old_mu = 0
            old_h = 4
            old_k = 0
            old_l = 0

        print("Enter primary-reflection angles:")
        delta = input("Delta = [{:6.2f}]: ".format(old_delta))
        if not delta:
            delta = old_delta
        th = input("Theta = [{:6.2f}]: ".format(old_th))
        if not th:
            th = old_th
        chi = input("Chi = [{:6.2f}]: ".format(old_chi))
        if not chi:
            chi = old_chi
        phi = input("Phi = [{:6.2f}]: ".format(old_phi))
        if not phi:
            phi = old_phi
        gamma = input("Gamma = [{:6.2f}]: ".format(old_gamma))
        if not gamma:
            gamma = old_gamma
        mu = input("Mu = [{:6.2f}]: ".format(old_mu))
        if not mu:
            mu = old_mu
        h = input("H = [{}]: ".format(old_h))
        if not h:
            h = old_h
        k = input("K = [{}]: ".format(old_k))
        if not k:
            k = old_k
        l = input("L = [{}]: ".format(old_l))
        if not l:
            l = old_l

    sample.add_reflection(
        float(h),
        float(k),
        float(l),
        position=_geom_.calc.Position(
            delta=float(delta),
            omega=float(th),
            chi=float(chi),
            phi=float(phi),
            gamma=float(gamma),
            mu=float(mu),
        ),
    )
    if len(orienting_refl) > 1:
        sample._orientation_reflections.pop(0)
    sample._orientation_reflections.insert(
        0, sample._sample.reflections_get()[-1]
    )

    if len(orienting_refl) > 1:
        print("Computing UB!")
        sample.compute_UB(
            sample._orientation_reflections[0],
            sample._orientation_reflections[1],
        )
        _geom_.forward(1, 0, 0)


def setor1(
    delta=None,
    th=None,
    chi=None,
    phi=None,
    gamma=None,
    mu=None,
    h=None,
    k=None,
    l=None,
):
    """
    Sets the primary secondary in hklpy.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    delta, th, chi, phi, gamma, mu : float, optional
        Values of motor positions for current reflection. If None, it will ask
        for it.
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    sample = _geom_.calc._sample
    ref = sample._sample.reflections_get()
    orienting_refl = sample._orientation_reflections
    if (
        not delta
        and not th
        and not chi
        and not phi
        and not gamma
        and not mu
        and not h
        and not k
        and not l
    ):
        if len(orienting_refl) > 1:
            for ref in sample._sample.reflections_get():
                if ref == orienting_refl[1]:
                    pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                    old_delta = pos[4]
                    old_th = pos[1]
                    old_chi = pos[2]
                    old_phi = pos[3]
                    old_gamma = pos[5]
                    old_mu = pos[0]
                    old_h, old_k, old_l = ref.hkl_get()

        else:
            old_delta = 60
            old_th = 30
            old_chi = 90
            old_phi = 0
            old_gamma = 0
            old_mu = 0
            old_h = 0
            old_k = 4
            old_l = 0

        print("Enter secondary-reflection angles:")
        delta = input("Delta = [{:6.2f}]: ".format(old_delta))
        if not delta:
            delta = old_delta
        th = input("Theta = [{:6.2f}]: ".format(old_th))
        if not th:
            th = old_th
        chi = input("Chi = [{:6.2f}]: ".format(old_chi))
        if not chi:
            chi = old_chi
        phi = input("Phi = [{:6.2f}]: ".format(old_phi))
        if not phi:
            phi = old_phi
        gamma = input("Gamma = [{:6.2f}]: ".format(old_gamma))
        if not gamma:
            gamma = old_gamma
        mu = input("Mu = [{:6.2f}]: ".format(old_mu))
        if not mu:
            mu = old_mu
        h = input("H = [{}]: ".format(old_h))
        if not h:
            h = old_h
        k = input("K = [{}]: ".format(old_k))
        if not k:
            k = old_k
        l = input("L = [{}]: ".format(old_l))
        if not l:
            l = old_l

    sample.add_reflection(
        float(h),
        float(k),
        float(l),
        position=_geom_.calc.Position(
            delta=float(delta),
            omega=float(th),
            chi=float(chi),
            phi=float(phi),
            gamma=float(gamma),
            mu=float(mu),
        ),
    )
    if len(orienting_refl) > 1:
        sample._orientation_reflections.pop(1)
    sample._orientation_reflections.insert(
        1, sample._sample.reflections_get()[-1]
    )

    if len(orienting_refl) > 1:
        print("Computing UB!")
        sample.compute_UB(
            sample._orientation_reflections[0],
            sample._orientation_reflections[1],
        )
        _geom_.forward(1, 0, 0)


def set_orienting():
    """
    Change the primary secondary orienting reflections to existing reflecitons
    in reflection list in hklpy.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.
    """
    sample = _geom_.calc._sample
    ref = sample._sample.reflections_get()
    orienting_refl = sample._orientation_reflections
    print(
        "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
            "#",
            "H",
            "K",
            "L",
            "Delta",
            "Theta",
            "Chi",
            "Phi",
            "Gamma",
            "Mu",
            "orienting",
        )
    )
    for i, ref in enumerate(sample._sample.reflections_get()):
        if orienting_refl[0] == ref:
            or0_old = i
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            print(
                "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                    i,
                    int(h),
                    int(k),
                    int(l),
                    pos[4],
                    pos[1],
                    pos[2],
                    pos[3],
                    pos[5],
                    pos[0],
                    "first",
                )
            )
        elif orienting_refl[1] == ref:
            or1_old = i
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            print(
                "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                    i,
                    int(h),
                    int(k),
                    int(l),
                    pos[4],
                    pos[1],
                    pos[2],
                    pos[3],
                    pos[5],
                    pos[0],
                    "second",
                )
            )
        else:
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            print(
                "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}  ".format(
                    i,
                    int(h),
                    int(k),
                    int(l),
                    pos[4],
                    pos[1],
                    pos[2],
                    pos[3],
                    pos[5],
                    pos[0],
                )
            )

    or0 = input("\nFirst orienting ({})? ".format(or0_old))
    if not or0:
        or0 = or0_old
    or1 = input("Second orienting ({})? ".format(or1_old))
    if not or1:
        or1 = or1_old
    sample._orientation_reflections.pop(0)
    sample._orientation_reflections.insert(
        0, sample._sample.reflections_get()[int(or0)]
    )
    sample._orientation_reflections.pop(1)
    sample._orientation_reflections.insert(
        1, sample._sample.reflections_get()[int(or1)]
    )
    print("Computing UB!")
    sample.compute_UB(
        sample._orientation_reflections[0], sample._orientation_reflections[1]
    )
    _geom_.forward(1, 0, 0)


def list_orienting(all_samples=False):
    """
    Prints the two reflections used in the UB matrix.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    all_samples : boolean, optional
        If True, it will print the reflections of all samples, if False, only of
        the current one.
    """
    _check_geom_selected()
    if all_samples:
        samples = _geom_.calc._samples.values()
    else:
        samples = [_geom_.calc._sample]
    for sample in samples:
        orienting_refl = sample._orientation_reflections
        print(
            "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
                "#",
                "H",
                "K",
                "L",
                "Delta",
                "Theta",
                "Chi",
                "Phi",
                "Gamma",
                "Mu",
                "orienting",
            )
        )
        for i, ref in enumerate(sample._sample.reflections_get()):
            if orienting_refl[0] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                        "first",
                    )
                )
            elif orienting_refl[1] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[4],
                        pos[1],
                        pos[2],
                        pos[3],
                        pos[5],
                        pos[0],
                        "second",
                    )
                )


def or0(h=None, k=None, l=None):
    """
    Sets the primary orientation in hklpy using the current motor positions.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections
    if not h and not k and not l:
        if len(orienting_refl) > 1:
            for ref in sample._sample.reflections_get():
                if ref == orienting_refl[0]:
                    hr, kr, lr = ref.hkl_get()
        else:
            hr = 2
            kr = 0
            lr = 0

        h = input("H ({})? ".format(hr)) if not h else h
        if not h:
            h = hr
        k = input("K ({})? ".format(kr)) if not k else k
        if not k:
            k = kr
        l = input("L ({})? ".format(lr)) if not l else l
        if not l:
            l = lr

    sample.add_reflection(
        float(h),
        float(k),
        float(l),
        position=_geom_.calc.Position(
            delta=_geom_.delta.get().user_readback,
            omega=_geom_.omega.get().user_readback,
            chi=_geom_.chi.get().user_readback,
            phi=_geom_.phi.get().user_readback,
            gamma=_geom_.gamma.get().user_readback,
            mu=_geom_.mu.get().user_readback,
        ),
    )

    if len(orienting_refl) > 1:
        sample._orientation_reflections.pop(0)
    sample._orientation_reflections.insert(
        0, sample._sample.reflections_get()[-1]
    )

    if len(orienting_refl) > 1:
        print("Computing UB!")
        sample.compute_UB(
            sample._orientation_reflections[0],
            sample._orientation_reflections[1],
        )
        _geom_.forward(1, 0, 0)


def or1(h=None, k=None, l=None):
    """
    Sets the secondary orientation in hklpy using the current motor positions.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections
    if not h and not k and not l:
        if len(orienting_refl) > 1:
            for ref in sample._sample.reflections_get():
                if ref == orienting_refl[1]:
                    hr, kr, lr = ref.hkl_get()
        else:
            hr = 0
            kr = 2
            lr = 0

        h = input("H ({})? ".format(hr)) if not h else h
        if not h:
            h = hr
        k = input("K ({})? ".format(kr)) if not k else k
        if not k:
            k = kr
        l = input("L ({})? ".format(lr)) if not l else l
        if not l:
            l = lr

    sample.add_reflection(
        float(h),
        float(k),
        float(l),
        position=_geom_.calc.Position(
            delta=_geom_.delta.get().user_readback,
            omega=_geom_.omega.get().user_readback,
            chi=_geom_.chi.get().user_readback,
            phi=_geom_.phi.get().user_readback,
            gamma=_geom_.gamma.get().user_readback,
            mu=_geom_.mu.get().user_readback,
        ),
    )

    if len(orienting_refl) > 1:
        sample._orientation_reflections.pop(1)
    sample._orientation_reflections.insert(
        1, sample._sample.reflections_get()[-1]
    )

    if len(orienting_refl) > 1:
        print("Computing UB!")
        sample.compute_UB(
            sample._orientation_reflections[0],
            sample._orientation_reflections[1],
        )
        _geom_.forward(1, 0, 0)


def compute_UB():
    """
    Calculates the UB matrix.

    This fixes one issue with the hklpy calc_UB in that using wh() right after
    will not work, it needs to run one calculation first.

    Parameters
    ----------
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    sample = _geom_.calc._sample
    print("Computing UB!")
    calc_UB(
        sample._orientation_reflections[0], sample._orientation_reflections[1]
    )
    _geom_.forward(1, 0, 0)


def calc_UB(r1, r2, wavelength=None, output=False):
    """
    Compute the UB matrix with two reflections.

    Parameters
    ----------
    r1, r2 : hklpy reflections
        Orienting reflections from hklpy.
    wavelength : float, optional
        This is not used...
    output : boolean
        Toggle to decide whether to print the UB matrix.
    """
    _check_geom_selected()
    _geom_.calc.sample.compute_UB(r1, r2)
    if output:
        print(_geom_.calc.sample.UB)


def setmode(mode=None):
    """
    Set the mode of the currently selected diffractometer.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    mode : string, optional
        Mode to be selected. If None, it will ask.
    """
    current_mode = _geom_.calc.engine.mode
    for index, item in enumerate(_geom_.calc.engine.modes):
        print("{:2d}. {}".format(index + 1, item))
        if current_mode == item:
            current_index = index
    if mode:
        _geom_.calc.engine.mode = _geom_.calc.engine.modes[int(mode) - 1]
        print("\nSet mode to {}".format(mode))
    else:
        mode = input("\nMode ({})? ".format(current_index + 1))
        if not mode:
            mode = current_index - 1
        _geom_.calc.engine.mode = _geom_.calc.engine.modes[int(mode) - 1]


def ca(h, k, l):
    """
    Calculate the motors position of a reflection.

    Parameters
    ----------
    h, k, l : float
        H, K, and L values.
    """

    pos = cahkl(h, k, l)
    print("\n   Calculated Positions:")
    print(
        "\n   H K L = {:5f} {:5f} {:5f}".format(
            h,
            k,
            l,
        )
    )
    print(
        "\n   Lambda (Energy) = {:6.4f} \u212B ({:6.4f} keV)".format(
            _geom_.calc.wavelength, _geom_.calc.energy
        )
    )
    print(
        "\n   {:8}{:8}{:8}{:8}{:8}{:8}".format(
            "Delta", "Theta", "Chi", "Phi", "Gamma", "Mu"
        )
    )
    print(
        "{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}".format(
            pos[4],
            pos[1],
            pos[2],
            pos[3],
            pos[5],
            pos[0],
        )
    )


def br(h, k, l):
    """
    Move the motors to a reciprocal space point.

    Parameters
    ----------
    h, k, l : float
        H, K, and L values.

    Returns
    -------
    Generator for the bluesky Run Engine.
    """
    yield from bps.mv(
        _geom_.h, float(h), _geom_.k, float(k), _geom_.l, float(l)
    )


def uan(delta=None, th=None):
    """
    Moves the delta and theta motors.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    delta, th: float, optional??
        Delta and th motor angles to be moved to.

    Returns
    -------
    Generator for the bluesky Run Engine.
    """
    if not delta or not th:
        raise ValueError("Usage: uan(delta,th)")
    else:
        print("Moving to (delta,th)=({},{})".format(delta, th))
        yield from bps.mv(_geom_.delta, delta, _geom_.omega, th)


def wh():
    """
    Retrieve information on the current reciprocal space position.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.
    """
    print(
        "\n   H K L = {:5f} {:5f} {:5f}".format(
            _geom_.calc.engine.pseudo_axes["h"],
            _geom_.calc.engine.pseudo_axes["k"],
            _geom_.calc.engine.pseudo_axes["l"],
        )
    )
    print(
        "\n   Lambda (Energy) = {:6.4f} \u212B ({:6.4f} keV)".format(
            _geom_.calc.wavelength, _geom_.calc.energy
        )
    )
    print(
        "\n   {:8}{:8}{:8}{:8}{:8}{:8}".format(
            "Delta", "Theta", "Chi", "Phi", "Gamma", "Mu"
        )
    )
    print(
        "{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}".format(
            _geom_.delta.get()[0],
            _geom_.omega.get()[0],
            _geom_.chi.get()[0],
            _geom_.phi.get()[0],
            _geom_.gamma.get()[0],
            _geom_.mu.get()[0],
        )
    )


def setlat(a=None, b=None, c=None, alpha=None, beta=None, gamma=None):
    """
    Set the lattice constants.

    Parameters
    ----------
    a, b, c, alpha, beta, gamma : float, optional
        Lattice constants. If None, it will ask for input.
    """

    current_sample = _geom_.calc.sample_name
    sample = _geom_.calc._samples[current_sample]
    lattice = [getattr(sample.lattice, parm) for parm in sample.lattice._fields]

    a = input("Lattice a ({})? ".format(lattice[0])) if not a else a
    if not a:
        a = lattice[0]
    b = input("Lattice b ({})? ".format(lattice[1])) if not b else b
    if not b:
        b = lattice[1]
    c = input("Lattice c ({})? ".format(lattice[2])) if not c else c
    if not c:
        c = lattice[2]
    alpha = (
        input("Lattice alpha ({})? ".format(lattice[3])) if not alpha else alpha
    )
    if not alpha:
        alpha = lattice[3]
    beta = input("Lattice beta ({})? ".format(lattice[4])) if not beta else beta
    if not beta:
        beta = lattice[4]
    gamma = (
        input("Lattice gamma ({})? ".format(lattice[5])) if not gamma else gamma
    )
    if not gamma:
        gamma = lattice[5]

    _geom_.calc.sample.lattice = (
        float(a),
        float(b),
        float(c),
        float(alpha),
        float(beta),
        float(gamma),
    )
    orienting_refl = sample._orientation_reflections
    if len(orienting_refl) > 1:
        print("Compute UB!")
        sample.compute_UB(
            sample._orientation_reflections[0],
            sample._orientation_reflections[1],
        )
        _geom_.forward(1, 0, 0)
