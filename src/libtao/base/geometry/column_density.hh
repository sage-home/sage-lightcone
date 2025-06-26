#ifndef tao_modules_columnDensity_hh
#define tao_modules_columnDensity_hh

#include <libhpc/numerics/simpson.hh>

namespace tao {
    class column_density {
      public:
        static real_type calculate_weighting_integral(real_type smoothing,
                                                      real_type distance_from_beam,
                                                      real_type distance_from_origin) {
            // Because the simulation is 3d, this will *always* be 3.
            // Apparently some simulations are done in 2D where this would
            // change.
            int dimensions = 3;
            if(distance_from_beam > smoothing) {
                real_type result = 0;
                return result;
            }
            return intKernel(distance_from_beam, smoothing, dimensions);
        }

      protected:
        static real_type get_weight(real_type normalized_impact_point, int dim) {
            // Get the dimensionless value of the kernel in units of smoothing length.
            // Springel 2005 Gadget2 paper Eqn 4 weighting function
            real_type f;
            if(normalized_impact_point < 0.5) {
                f = 1.0 - 6.0 * pow(normalized_impact_point, 2) + 6 * pow(normalized_impact_point, 3);
            } else if(normalized_impact_point < 1) {
                f = 2.0 * pow(1.0 - normalized_impact_point, 3);
            } else {
                f = 0.0;
            }

            // To get agreement with Price 2005 and Springel 2005 which has half
            // the smoothing length definition
            f *= pow(2, dim);

            // From Price 2005 this is the normalisation constants given in Eqn 5
            real_type sigma;
            switch(dim) {
                case 1:
                    sigma = 2.0 / 2.3;
                    break;
                case 2:
                    sigma = 10.0 / (7 * M_PI);
                    break;
                case 3:
                    sigma = 1.0 / M_PI;
                    break;
                default:
                    throw "Invalid dimensions";
            }

            return f * sigma;
        }

        struct get_weight_func {
            // To get a wrapped weighting function with constant impact point
            // and constant dimensions, but a varying distance into the particle along
            // the ray / pencilbeam.
            typedef real_type value_type;

            get_weight_func(real_type impact_point, int dim) : impact_point(impact_point), dim(dim) {
            }

            real_type operator()(real_type distance) const {
                return get_weight(sqrt(pow(distance, 2) + pow(impact_point, 2)), dim);
            }

            real_type impact_point;
            int       dim;
        };

        static real_type intKernel(real_type impact_point, real_type smoothing_length, int dim, int kernelspline = 1) {
            // Price's eqn 30
            // Integrate the Kernel along ray axis at minimum distance (impact_point)
            // from origin

            // Put in units of smoothing length
            impact_point /= smoothing_length;

            // Limit is from impact point to edge of particle, along axis at right angles to
            // the particle's radial line to the impact point.
            // As we've normalized our impact_point by smoothing,
            // hypotenuse of the triangle (radius of particle) is 1 (the kernelspline).
            real_type zmax_2D = sqrt(pow(kernelspline, 2) - pow(impact_point, 2));
            // Integral is symmetric about zero,
            // so (-zmax_2D to 0) + (0 to zmax_2D) is just double (0 to zmax_2D)
            real_type integrated_los = 2.0 * hpc::num::simpson(
                                                 // Get the weighting function as a function of distance
                                                 get_weight_func(impact_point, dim),
                                                 // Limits
                                                 0,
                                                 zmax_2D,
                                                 // TODO Number of points to split into - Accuracy
                                                 50);

            switch(dim) {
                case 3:
                    // 3D particle, integrated along line of sight means it returns
                    // a 2D value, Prince 2005 Eqn 29
                    integrated_los /= pow(smoothing_length, 2);
                    break;
                case 2:
                    // 2D particle (perhaps a test accretion disk simulation?)
                    // that wants integrated LoS so now 1D
                    integrated_los /= smoothing_length;
                    break;
                default:
                    // Can't integrate the 1D test simulation along the line of sight
                    throw "Invalid dimensions";
            }
            return integrated_los;
        }
    };
} // namespace tao

#endif
