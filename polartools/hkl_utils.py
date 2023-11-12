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
    ~del_reflection
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
    ~update_lattice
    ~read_config
    ~write_config
"""

import bluesky.plan_stubs as bps
import pathlib

try:
    from hkl import cahkl
    from hkl.user import _check_geom_selected, current_diffractometer
    from hkl.configuration import DiffractometerConfiguration
except ModuleNotFoundError:
    print("Some module did not load!")

path_startup = pathlib.Path("startup_experiment.py")


def sampleChange(sample_key=None):
    """
    Change selected sample in hklpy.

    Parameters
    ----------
    sample_key : string, optional
        Name of the sample as set in hklpy. If None it will ask for which
        sample.
    """
    _geom_ = current_diffractometer()

    if sample_key is None:
        d = _geom_.calc._samples.keys()
        print("Sample keys:", list(d))
        sample_key = (
            input("\nEnter sample key [{}]: ".format(_geom_.calc.sample.name))
            or _geom_.calc.sample.name
        )
    try:
        _geom_.calc.sample = _geom_.calc._samples[
            sample_key
        ]  # define the current sample
        print("\nCurrent sample: " + _geom_.calc.sample.name)
        # to be done: check if orienting reflections exist
        compute_UB()

    except KeyError:
        print("Not a valid sample key")


def sampleList():
    """List all samples currently defined in hklpy; specify  current one."""

    _geom_ = current_diffractometer()

    samples = _geom_.calc._samples
    for x in list(samples.keys())[1:]:
        orienting_refl = samples[x]._orientation_reflections
        print("\nSample = {}".format(x))
        print("Lattice:", end=" ")
        print(*samples[x].lattice._fields, sep=", ", end=" = ")
        print(*samples[x].lattice, sep=", ")
        for ref in samples[x]._sample.reflections_get():
            if orienting_refl[0] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                if len(_geom_.calc.physical_axes) == 6:
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
                elif len(_geom_.calc.physical_axes) == 4:
                    print(
                        "{:>3}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}   ".format(
                            "or0",
                            int(h),
                            int(k),
                            int(l),
                            pos[3],
                            pos[0],
                            pos[1],
                            pos[2],
                        )
                    )
                else:
                    raise ValueError(
                        "Geometry {} not supported.".format(_geom_.name)
                    )

            elif orienting_refl[1] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                if len(_geom_.calc.physical_axes) == 6:
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
                elif len(_geom_.calc.physical_axes) == 4:
                    print(
                        "{:>3}{:>4}{:>3}{:>3}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}  ".format(
                            "or1",
                            int(h),
                            int(k),
                            int(l),
                            pos[3],
                            pos[0],
                            pos[1],
                            pos[2],
                        )
                    )
                else:
                    raise ValueError(
                        "Geometry {} not supported.".format(_geom_.name)
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
    _geom_ = current_diffractometer()
    _check_geom_selected()

    if all_samples:
        samples = _geom_.calc._samples.values()
    else:
        samples = [_geom_.calc._sample]
    for sample in samples:
        print("Sample: {}".format(sample.name))
        orienting_refl = sample._orientation_reflections
        if len(_geom_.calc.physical_axes) == 6:
            print(
                "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
                    "#",
                    "H",
                    "K",
                    "L",
                    "Delta",
                    "Eta",
                    "Chi",
                    "Phi",
                    "Nu",
                    "Mu",
                    "orienting",
                )
            )
        elif len(_geom_.calc.physical_axes) == 4:
            print(
                "\n{:>2}{:>4}{:>3}{:>3}{:>12}{:>9}{:>9}{:>9}   {:<12}".format(
                    "#",
                    "H",
                    "K",
                    "L",
                    "Two Theta",
                    "Theta",
                    "Chi",
                    "Phi",
                    "orienting",
                )
            )
        else:
            raise ValueError("Geometry {} not supported.".format(_geom_.name))

        for i, ref in enumerate(sample._sample.reflections_get()):
            if orienting_refl[0] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                if len(_geom_.calc.physical_axes) == 6:
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
                elif len(_geom_.calc.physical_axes) == 4:
                    print(
                        "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                            i,
                            int(h),
                            int(k),
                            int(l),
                            pos[3],
                            pos[0],
                            pos[1],
                            pos[2],
                            "first",
                        )
                    )
                else:
                    raise ValueError(
                        "Geometry {} not supported.".format(_geom_.name)
                    )
            elif orienting_refl[1] == ref:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                if len(_geom_.calc.physical_axes) == 6:
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
                elif len(_geom_.calc.physical_axes) == 4:
                    print(
                        "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                            i,
                            int(h),
                            int(k),
                            int(l),
                            pos[3],
                            pos[0],
                            pos[1],
                            pos[2],
                            "second",
                        )
                    )
                else:
                    raise ValueError(
                        "Geometry {} not supported.".format(_geom_.name)
                    )
            else:
                h, k, l = ref.hkl_get()
                pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                if len(_geom_.calc.physical_axes) == 6:
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
                elif len(_geom_.calc.physical_axes) == 4:
                    print(
                        "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f} ".format(
                            i,
                            int(h),
                            int(k),
                            int(l),
                            pos[3],
                            pos[0],
                            pos[1],
                            pos[2],
                        )
                    )
                else:
                    raise ValueError(
                        "Geometry {} not supported.".format(_geom_.name)
                    )
        if len(samples) > 1 and all_samples:
            print(
                "============================================================================"
            )


def or_swap():
    """Swaps the two orientation reflections in hklpy."""
    _geom_ = current_diffractometer()
    sample = _geom_.calc._sample
    sample.swap_orientation_reflections()
    list_reflections()
    print("Computing UB!")
    sample.compute_UB(
        sample._orientation_reflections[0], sample._orientation_reflections[1]
    )
    # calculate reflection to update values. solutions for 1st orienting reflections always exist
    first_orienting = sample._orientation_reflections[0].hkl_get()
    _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


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
    _geom_ = current_diffractometer()
    _check_geom_selected()
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections

    if len(_geom_.calc.physical_axes) == 6 and len(args) == 9:
        delta, th, chi, phi, gamma, mu, h, k, l = args
    elif len(_geom_.calc.physical_axes) == 4 and len(args) == 7:
        delta, th, chi, phi, h, k, l = args
    else:
        if len(orienting_refl) > 1:
            for ref in sample._sample.reflections_get():
                if (
                    ref == orienting_refl[0]
                    and len(_geom_.calc.physical_axes) == 6
                ):
                    pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                    old_delta = pos[4]
                    old_th = pos[1]
                    old_chi = pos[2]
                    old_phi = pos[3]
                    old_gamma = pos[5]
                    old_mu = pos[0]
                    old_h, old_k, old_l = ref.hkl_get()
                elif (
                    ref == orienting_refl[0]
                    and len(_geom_.calc.physical_axes) == 4
                ):
                    pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                    old_delta = pos[3]
                    old_th = pos[0]
                    old_chi = pos[1]
                    old_phi = pos[2]
                    old_h, old_k, old_l = ref.hkl_get()

        else:
            old_delta = 60
            old_th = 30
            old_chi = 90
            old_phi = 0
            old_h = 4
            old_k = 0
            old_l = 0
            old_gamma = 0
            old_mu = 0

        print("Enter primary-reflection angles:")
        delta = input("Delta = [{:6.2f}]: ".format(old_delta)) or old_delta
        th = input("Eta = [{:6.2f}]: ".format(old_th)) or old_th
        chi = input("Chi = [{:6.2f}]: ".format(old_chi)) or old_chi
        phi = input("Phi = [{:6.2f}]: ".format(old_phi)) or old_phi
        if len(_geom_.calc.physical_axes) == 6:
            gamma = input("Nu = [{:6.2f}]: ".format(old_gamma)) or old_gamma
            mu = input("Mu = [{:6.2f}]: ".format(old_mu)) or old_mu
        h = input("H = [{}]: ".format(old_h)) or old_h
        k = input("K = [{}]: ".format(old_k)) or old_k
        l = input("L = [{}]: ".format(old_l)) or old_l

    if len(_geom_.calc.physical_axes) == 6:
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
    elif len(_geom_.calc.physical_axes) == 4:
        sample.add_reflection(
            float(h),
            float(k),
            float(l),
            position=_geom_.calc.Position(
                tth=float(delta),
                omega=float(th),
                chi=float(chi),
                phi=float(phi),
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
        # calculate reflection to update values. solutions for 1st orienting reflections always exist
        first_orienting = sample._orientation_reflections[0].hkl_get()
        _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


def setor1(*args):
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

    _geom_ = current_diffractometer()
    _check_geom_selected()
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections

    if len(_geom_.calc.physical_axes) == 6 and len(args) == 9:
        delta, th, chi, phi, gamma, mu, h, k, l = args
    elif len(_geom_.calc.physical_axes) == 4 and len(args) == 7:
        delta, th, chi, phi, h, k, l = args
    else:
        if len(orienting_refl) > 1:
            for ref in sample._sample.reflections_get():
                if (
                    ref == orienting_refl[1]
                    and len(_geom_.calc.physical_axes) == 6
                ):
                    pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                    old_delta = pos[4]
                    old_th = pos[1]
                    old_chi = pos[2]
                    old_phi = pos[3]
                    old_gamma = pos[5]
                    old_mu = pos[0]
                    old_h, old_k, old_l = ref.hkl_get()
                elif (
                    ref == orienting_refl[1]
                    and len(_geom_.calc.physical_axes) == 4
                ):
                    pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
                    old_delta = pos[3]
                    old_th = pos[0]
                    old_chi = pos[1]
                    old_phi = pos[2]
                    old_h, old_k, old_l = ref.hkl_get()

        else:
            old_delta = 60
            old_th = 30
            old_chi = 0
            old_phi = 0
            old_h = 0
            old_k = 4
            old_l = 0
            old_gamma = 0
            old_mu = 0

        print("Enter secondary-reflection angles:")
        delta = input("Delta = [{:6.2f}]: ".format(old_delta)) or old_delta
        th = input("Eta = [{:6.2f}]: ".format(old_th)) or old_th
        chi = input("Chi = [{:6.2f}]: ".format(old_chi)) or old_chi
        phi = input("Phi = [{:6.2f}]: ".format(old_phi)) or old_phi
        if len(_geom_.calc.physical_axes) == 6:
            gamma = input("Nu = [{:6.2f}]: ".format(old_gamma)) or old_gamma
            mu = input("Mu = [{:6.2f}]: ".format(old_mu)) or old_mu
        h = input("H = [{}]: ".format(old_h)) or old_h
        k = input("K = [{}]: ".format(old_k)) or old_k
        l = input("L = [{}]: ".format(old_l)) or old_l

    if len(_geom_.calc.physical_axes) == 6:
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
    elif len(_geom_.calc.physical_axes) == 4:
        sample.add_reflection(
            float(h),
            float(k),
            float(l),
            position=_geom_.calc.Position(
                tth=float(delta),
                omega=float(th),
                chi=float(chi),
                phi=float(phi),
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
        # calculate reflection to update values. solutions for 1st orienting reflections always exist
        first_orienting = sample._orientation_reflections[0].hkl_get()
        _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


def set_orienting():
    """
    Change the primary secondary orienting reflections to existing reflecitons
    in reflection list in hklpy.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.
    """
    _geom_ = current_diffractometer()
    _check_geom_selected()
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections
    if len(_geom_.calc.physical_axes) == 6:
        print(
            "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
                "#",
                "H",
                "K",
                "L",
                "Delta",
                "Eta",
                "Chi",
                "Phi",
                "Nu",
                "Mu",
                "orienting",
            )
        )
    elif len(_geom_.calc.physical_axes) == 4:
        print(
            "\n{:>2}{:>4}{:>3}{:>3}{:>12}{:>9}{:>9}{:>9}   {:<12}".format(
                "#",
                "H",
                "K",
                "L",
                "Two Theta",
                "Theta",
                "Chi",
                "Phi",
                "orienting",
            )
        )
    else:
        raise ValueError("Geometry {} not supported.".format(_geom_.name))

    for i, ref in enumerate(sample._sample.reflections_get()):
        if orienting_refl[0] == ref:
            or0_old = i
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                        "first",
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )
        elif orienting_refl[1] == ref:
            or1_old = i
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                        "second",
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )
        else:
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )

    or0 = input("\nFirst orienting ({})? ".format(or0_old)) or or0_old
    or1 = input("Second orienting ({})? ".format(or1_old)) or or1_old
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
    # calculate reflection to update values. solutions for 1st orienting reflections always exist
    first_orienting = sample._orientation_reflections[0].hkl_get()
    _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


def del_reflection():
    """
    Delete existing reflection from in reflection list in hklpy.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.
    """
    _geom_ = current_diffractometer()
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections
    if len(_geom_.calc.physical_axes) == 6:
        print(
            "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
                "#",
                "H",
                "K",
                "L",
                "Delta",
                "Eta",
                "Chi",
                "Phi",
                "Nu",
                "Mu",
                "orienting",
            )
        )
    elif len(_geom_.calc.physical_axes) == 4:
        print(
            "\n{:>2}{:>4}{:>3}{:>3}{:>12}{:>9}{:>9}{:>9}   {:<12}".format(
                "#",
                "H",
                "K",
                "L",
                "Two Theta",
                "Theta",
                "Chi",
                "Phi",
                "orienting",
            )
        )
    else:
        raise ValueError("Geometry {} not supported.".format(_geom_.name))

    for i, ref in enumerate(sample._sample.reflections_get()):
        if orienting_refl[0] == ref:
            or0_old = i
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                        "first",
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )
        elif orienting_refl[1] == ref:
            or1_old = i
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                        "second",
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )
        else:
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )

    remove = input("\nRemove reflection # ")
    if not remove:
        print("No reflection removed")
    elif int(remove) == or0_old or int(remove) == or1_old:
        print("Orienting reflection not removable!")
        print(
            "Use 'set_orienting()' first to select different orienting reflection."
        )
    else:
        sample._sample.del_reflection(
            sample._sample.reflections_get()[int(remove)]
        )


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
    _geom_ = current_diffractometer()
    _check_geom_selected()
    if all_samples:
        samples = _geom_.calc._samples.values()
    else:
        samples = [_geom_.calc._sample]
    for sample in samples:
        orienting_refl = sample._orientation_reflections
        if len(_geom_.calc.physical_axes) == 6:
            print(
                "\n{:>2}{:>4}{:>3}{:>3}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}   {:<12}".format(
                    "#",
                    "H",
                    "K",
                    "L",
                    "Delta",
                    "Eta",
                    "Chi",
                    "Phi",
                    "Nu",
                    "Mu",
                    "orienting",
                )
            )
        elif len(_geom_.calc.physical_axes) == 4:
            print(
                "\n{:>2}{:>4}{:>3}{:>3}{:>12}{:>9}{:>9}{:>9}   {:<12}".format(
                    "#",
                    "H",
                    "K",
                    "L",
                    "Two Theta",
                    "Theta",
                    "Chi",
                    "Phi",
                    "orienting",
                )
            )
        else:
            raise ValueError("Geometry {} not supported.".format(_geom_.name))

    for i, ref in enumerate(sample._sample.reflections_get()):
        if orienting_refl[0] == ref:
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                        "first",
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
                )
        elif orienting_refl[1] == ref:
            h, k, l = ref.hkl_get()
            pos = ref.geometry_get().axis_values_get(_geom_.calc._units)
            if len(_geom_.calc.physical_axes) == 6:
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
            elif len(_geom_.calc.physical_axes) == 4:
                print(
                    "{:>2}{:>4}{:>3}{:>3}{:>12.3f}{:>9.3f}{:>9.3f}{:>9.3f}   {:<12} ".format(
                        i,
                        int(h),
                        int(k),
                        int(l),
                        pos[3],
                        pos[0],
                        pos[1],
                        pos[2],
                        "second",
                    )
                )
            else:
                raise ValueError(
                    "Geometry {} not supported.".format(_geom_.name)
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
    _geom_ = current_diffractometer()
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
        h = (input("H ({})? ".format(hr)) if not h else h) or hr
        k = (input("K ({})? ".format(kr)) if not k else k) or kr
        l = (input("L ({})? ".format(lr)) if not l else l) or lr
    if len(_geom_.calc.physical_axes) == 6:
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
    elif len(_geom_.calc.physical_axes) == 4:
        sample.add_reflection(
            float(h),
            float(k),
            float(l),
            position=_geom_.calc.Position(
                tth=_geom_.tth.get().user_readback,
                omega=_geom_.omega.get().user_readback,
                chi=_geom_.chi.get().user_readback,
                phi=_geom_.phi.get().user_readback,
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
        # calculate reflection to update values. solutions for 1st orienting reflections always exist
        first_orienting = sample._orientation_reflections[0].hkl_get()
        _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


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
    _geom_ = current_diffractometer()
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
        h = (input("H ({})? ".format(hr)) if not h else h) or hr
        k = (input("K ({})? ".format(kr)) if not k else k) or kr
        l = (input("L ({})? ".format(lr)) if not l else l) or lr
    if len(_geom_.calc.physical_axes) == 6:
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
    elif len(_geom_.calc.physical_axes) == 4:
        sample.add_reflection(
            float(h),
            float(k),
            float(l),
            position=_geom_.calc.Position(
                tth=_geom_.tth.get().user_readback,
                omega=_geom_.omega.get().user_readback,
                chi=_geom_.chi.get().user_readback,
                phi=_geom_.phi.get().user_readback,
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
        # calculate reflection to update values. solutions for 1st orienting reflections always exist
        first_orienting = sample._orientation_reflections[0].hkl_get()
        _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


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
    _geom_ = current_diffractometer()
    sample = _geom_.calc._sample
    print("Computing UB!")
    calc_UB(
        sample._orientation_reflections[0], sample._orientation_reflections[1]
    )
    # calculate reflection to update values. solutions for 1st orienting reflections always exist
    first_orienting = sample._orientation_reflections[0].hkl_get()
    _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


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
    _geom_ = current_diffractometer()
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
    _geom_ = current_diffractometer()
    current_mode = _geom_.calc.engine.mode
    for index, item in enumerate(_geom_.calc.engine.modes):
        print("{:2d}. {}".format(index + 1, item))
        if current_mode == item:
            current_index = index
    if mode:
        _geom_.calc.engine.mode = _geom_.calc.engine.modes[int(mode) - 1]
        print("\nSet mode to {}".format(mode))
    else:
        mode = input("\nMode ({})? ".format(current_index + 1)) or (
            current_index + 1
        )
        _geom_.calc.engine.mode = _geom_.calc.engine.modes[int(mode) - 1]


def ca(h, k, l):
    """
    Calculate the motors position of a reflection.

    Parameters
    ----------
    h, k, l : float
        H, K, and L values.
    """
    _geom_ = current_diffractometer()
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
    if len(_geom_.calc.physical_axes) == 6:
        print(
            "\n{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}".format(
                "Delta", "Eta", "Chi", "Phi", "Nu", "Mu"
            )
        )
        print(
            "{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}".format(
                pos[4],
                pos[1],
                pos[2],
                pos[3],
                pos[5],
                pos[0],
            )
        )
    elif len(_geom_.calc.physical_axes) == 4:
        print(
            "\n{:>9}{:>9}{:>9}{:>9}".format("Two Theta", "Theta", "Chi", "Phi")
        )
        print(
            "{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}".format(
                pos[3],
                pos[0],
                pos[1],
                pos[2],
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
    _geom_ = current_diffractometer()
    yield from bps.mv(
        _geom_.h, float(h), _geom_.k, float(k), _geom_.l, float(l)
    )
    pos = cahkl(h, k, l)
    if len(_geom_.calc.physical_axes) == 6:
        print(
            "\n Delta, Eta, Chi, Phi, Nu, Mu --> {:>6.3f}, {:>6.3f}, {:>6.3f}, {:>6.3f}, {:>6.3f}, {:>6.3f}".format(
                pos[4],
                pos[1],
                pos[2],
                pos[3],
                pos[5],
                pos[0],
            )
        )
    elif len(_geom_.calc.physical_axes) == 4:
        print(
            "Two Theta, Theta, Chi, Phi --> {:>6.3f}, {:>6.3f}, {:>6.3f}, {:>6.3f}".format(
                pos[3],
                pos[0],
                pos[1],
                pos[2],
            )
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
    _geom_ = current_diffractometer()
    if not delta or not th:
        raise ValueError("Usage: uan(delta/tth,th)")
    else:
        if len(_geom_.calc.physical_axes) == 6:
            print("Moving to (delta,eta)=({},{})".format(delta, th))
            yield from bps.mv(_geom_.delta, delta, _geom_.omega, th)
        elif len(_geom_.calc.physical_axes) == 4:
            print("Moving to (tth,th)=({},{})".format(delta, th))
            yield from bps.mv(_geom_.tth, delta, _geom_.omega, th)


def wh():
    """
    Retrieve information on the current reciprocal space position.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.
    """
    _geom_ = current_diffractometer()
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
    if len(_geom_.calc.physical_axes) == 6:
        print(
            "\n{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}".format(
                "Delta", "Eta", "Chi", "Phi", "Nu", "Mu"
            )
        )
        print(
            "{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}{:>9.3f}".format(
                _geom_.delta.get()[0],
                _geom_.omega.get()[0],
                _geom_.chi.get()[0],
                _geom_.phi.get()[0],
                _geom_.gamma.get()[0],
                _geom_.mu.get()[0],
            )
        )
    elif len(_geom_.calc.physical_axes) == 4:
        print(
            "\n{:>11}{:>9}{:>9}{:>9}".format("Two Theta", "Theta", "Chi", "Phi")
        )
        print(
            "{:>11.3f}{:>9.3f}{:>9.3f}{:>9.3f}".format(
                _geom_.tth.get()[0],
                _geom_.omega.get()[0],
                _geom_.chi.get()[0],
                _geom_.phi.get()[0],
            )
        )


def setlat(*args):
    """
    Set the lattice constants.

    Parameters
    ----------
    a, b, c, alpha, beta, gamma : float, optional
        Lattice constants. If None, it will ask for input.
    """

    _geom_ = current_diffractometer()
    current_sample = _geom_.calc.sample_name
    sample = _geom_.calc._samples[current_sample]
    lattice = [getattr(sample.lattice, parm) for parm in sample.lattice._fields]

    if len(args) == 6:
        a, b, c, alpha, beta, gamma = args
    elif len(args) == 0:
        a = (input("Lattice a ({})? ".format(lattice[0]))) or lattice[0]
        b = (input("Lattice b ({})? ".format(lattice[1]))) or lattice[1]
        c = (input("Lattice c ({})? ".format(lattice[2]))) or lattice[2]
        alpha = (input("Lattice alpha ({})? ".format(lattice[3]))) or lattice[3]
        beta = (input("Lattice beta ({})? ".format(lattice[4]))) or lattice[4]
        gamma = (input("Lattice gamma ({})? ".format(lattice[5]))) or lattice[5]
    else:
        raise ValueError(
            "either no arguments or a, b, c, alpha, beta, gamma need to be provided."
        )

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
        # calculate reflection to update values. solutions for 1st orienting reflections always exist
        first_orienting = sample._orientation_reflections[0].hkl_get()
        _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)


def update_lattice(lattice_constant=None):
    """
    Update lattice constants.

    Parameters
    ----------
    lattice_constant: string, optional
        a, b or c or auto (default)
    """

    _geom_ = current_diffractometer()
    current_sample = _geom_.calc.sample_name
    sample = _geom_.calc._samples[current_sample]
    lattice = [getattr(sample.lattice, parm) for parm in sample.lattice._fields]
    a = lattice[0]
    b = lattice[1]
    c = lattice[2]
    alpha = lattice[3]
    beta = lattice[4]
    gamma = lattice[5]
    hh = _geom_.calc.engine.pseudo_axes["h"]
    kk = _geom_.calc.engine.pseudo_axes["k"]
    ll = _geom_.calc.engine.pseudo_axes["l"]

    if round(hh) == 0 and round(kk) == 0:
        lattice_auto = "c"
    elif round(hh) == 0 and round(ll) == 0:
        lattice_auto = "b"
    elif round(kk) == 0 and round(ll) == 0:
        lattice_auto = "a"
    else:
        lattice_auto = None
    if not lattice_constant:
        lattice_constant = (
            input("Lattice parameter (a, b, or c or [auto])? ") or lattice_auto
        )
    else:
        print("Specify lattice parameter 'a', 'b' or 'c' or none")
    if lattice_constant == "a" and abs(round(hh)) > 0:
        a = a / hh * round(hh)
    elif lattice_constant == "b" and abs(round(kk)) > 0:
        b = b / kk * round(kk)
    elif lattice_constant == "c" and abs(round(ll)) > 0:
        c = c / ll * round(ll)
    else:
        raise ValueError("Auto calc not possible.")
    print("Refining lattice parameter {}".format(lattice_constant))
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
        # calculate reflection to update values. solutions for 1st orienting reflections always exist
        first_orienting = sample._orientation_reflections[0].hkl_get()
        _geom_.forward(first_orienting.h, first_orienting.k, first_orienting.l)
    print(
        "\n   H K L = {:5.4f} {:5.4f} {:5.4f}".format(
            _geom_.calc.engine.pseudo_axes["h"],
            _geom_.calc.engine.pseudo_axes["k"],
            _geom_.calc.engine.pseudo_axes["l"],
        )
    )
    lattice = [getattr(sample.lattice, parm) for parm in sample.lattice._fields]
    print(
        "   a, b, c, alpha, beta, gamma = {:5.4f} {:5.4f} {:5.4f} {:5.4f} {:5.4f} {:5.4f}".format(
            lattice[0],
            lattice[1],
            lattice[2],
            lattice[3],
            lattice[4],
            lattice[5],
        )
    )


def write_config(method="File", overwrite=False):
    """
    Write configuration from file in current directory.

    Parameters
    ----------
    method: string, optional
        right now only "File" possible, but later PV or other
    overwrite: Boolean, optional
        asks if existing file hould be overwritten
    """
    _geom_ = current_diffractometer()
    config = DiffractometerConfiguration(_geom_)
    if len(_geom_.calc.physical_axes) == 6:
        config_file = pathlib.Path("diffract-config.json")
    elif len(_geom_.calc.physical_axes) == 4:
        config_file = pathlib.Path("fourc-config.json")
    else:
        raise ValueError("Geometry not supported.")
    settings = config.export("json")
    if config_file.exists():
        if not overwrite:
            value = input("Overwrite existing configuration file (y/[n])? ")
            if value == "y":
                overwrite = True
        if overwrite:
            if method == "File":
                with open(config_file.name, "w") as f:
                    f.write(settings)
                print("Writing configuration file {}.".format(config_file.name))
    else:
        if method == "File":
            print("Writing configuration file.")
            with open(config_file.name, "w") as f:
                f.write(settings)


def read_config(method="File"):
    """
    Read configuration from file in current directory.

    Parameters
    ----------
    method: string, optional
        right now only "File" possible, but later PV or other
    """
    _geom_ = current_diffractometer()
    config = DiffractometerConfiguration(_geom_)
    sample = _geom_.calc._sample
    orienting_refl = sample._orientation_reflections

    if len(_geom_.calc.physical_axes) == 6:
        config_file = pathlib.Path("diffract-config.json")
    elif len(_geom_.calc.physical_axes) == 4:
        config_file = pathlib.Path("fourc-config.json")
    else:
        raise ValueError("Geometry not supported.")

    if config_file.exists():
        if method == "File":
            print("Read configuration file '{}'.".format(config_file.name))
            method = input("Method ([o]verwrite/[a]ppend)? ")
            if method == "a":
                config.restore(config_file, clear=False)
                print("Appending stored reflections.")
                if len(orienting_refl) > 1:
                    print("Compute UB!")
                    sample.compute_UB(
                        sample._orientation_reflections[0],
                        sample._orientation_reflections[1],
                    )
                first_orienting = sample._orientation_reflections[0].hkl_get()
                _geom_.forward(
                    first_orienting.h, first_orienting.k, first_orienting.l
                )
            elif method == "o":
                config.restore(config_file, clear=True)
                print("Overwriting all reflections with stored reflections.")
                if len(orienting_refl) > 1:
                    print("Compute UB!")
                    sample.compute_UB(
                        sample._orientation_reflections[0],
                        sample._orientation_reflections[1],
                    )
                first_orienting = sample._orientation_reflections[0].hkl_get()
                _geom_.forward(
                    first_orienting.h, first_orienting.k, first_orienting.l
                )
            else:
                raise ValueError("Either overwrite or append.")
        else:
            raise ValueError("Method not supported yet.")
    else:
        raise ValueError("No config file found.")
