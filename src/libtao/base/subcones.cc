#include "subcones.hh"
#include <libhpc/algorithm/ridders.hh>
#include <math.h>

namespace tao {

boost::optional<double> calc_subcone_angle(tao::lightcone const &lc) {
  EXCEPT(lc.max_ra() - lc.min_ra() <= M_PI / 2, "Subcone RA above 90 degrees.");
  EXCEPT(lc.max_dec() - lc.min_dec() <= M_PI / 2,
         "Subcone DEC above 90 degrees.");

  tao::simulation const *sim = lc.simulation();
  double b = sim->box_size();
  double ra = lc.max_ra() - lc.min_ra();
  double d0 = lc.min_dist();
  double d1 = lc.max_dist();

  // If the cone fits in one box then return zero.
  if (d1 - d0 * cos(ra) <= b)
    return 0.0 - lc.min_ra();

  // Use Ridder's method to find the optimal angle for
  // unique cones.
  auto res = hpc::ridders(
      [ra, d0, d1, b](double x) {
        double phi = ra + x;
        return b - d1 * (cos(x) - sin(x) / tan(phi));
      },
      0.5 * M_PI, 0.0);
  if (res != std::numeric_limits<double>::max())
    return res - lc.min_ra();
  else
    return boost::none;
}

unsigned calc_max_subcones(tao::lightcone const &lc) {
  auto theta = calc_subcone_angle(lc);
  if (!theta)
    return 0;
  double b = lc.simulation()->box_size();
  double d0 = lc.min_dist();
  double d1 = lc.max_dist();
  double phi = *theta + (lc.max_ra() - lc.min_ra());
  double h = d1 * sin(phi) - d0 * sin(*theta);
  double h_dec = d1 * sin((lc.max_dec() - lc.min_dec()));
  return (unsigned)(floor(b / h) * floor(b / h_dec));
}

} // namespace tao
