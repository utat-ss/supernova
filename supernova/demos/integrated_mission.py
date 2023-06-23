from supernova.api import propagate_orbit
from supernova.plotter import plot_from_array, plot_3d_from_array
from celest.satellite import Satellite
from celest.time import Time
from celest.coordinates import Coordinate, GroundLocation, GCRS, ITRS
from celest.encounter import generate_vtws
from celest import units as u

from datetime import datetime


def jd2000_to_datetime(jd2000: u.Quantity) -> datetime:
    return Time(jd2000).datetime


if __name__ == "__main__":
    days = 50

    t_span = [0, 86400 * days]

    y0 = [
        -4749231.102296294,
        -4975106.82469687,
        0.0,
        -719.6538503589323,
        686.980716081442,
        7561.282735263496,
    ]

    t, y = propagate_orbit("RK810", "simplified", t_span, y0, 1e-6)

    print(f"Steps taken: {len(t)}")
    # plot_from_array(t, y)
    # plot_3d_from_array(t, y)

    # Celest stuff
    julian = Time(t / 86400, offset=2414900)

    position = GCRS(julian.julian.data, y[:, 0], y[:, 1], y[:, 2], u.m)
    velocity = GCRS(julian.julian.data, y[:, 3], y[:, 4], y[:, 5], u.m / u.s)

    toronto = GroundLocation(
        latitude=43.6532,
        longitude=-79.3832,
        height=76,
        angular_unit=u.deg,
        length_unit=u.m,
    )

    # pos_aa = Coordinate(position).convert_to("AzEl", toronto)
    # pos_itrs = Coordinate(position).convert_to(ITRS)

    # plot_from_array(t, pos_itrs.data)

    satellite = Satellite(position=position, velocity=velocity)

    # # Generate ground location windows.
    downlinking_windows = generate_vtws(
        satellite=satellite, location=toronto, vis_threshold=10
    )

    for window in downlinking_windows:
        rise_time = window.rise_time

        print(jd2000_to_datetime(rise_time).strftime("%Y-%m-%d %H:%M:%S"))

        set_time = window.set_time

        print(jd2000_to_datetime(set_time).strftime("%Y-%m-%d %H:%M:%S"))
