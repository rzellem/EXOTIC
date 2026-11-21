// Claret radial quadrature, matching PyLightcurve 4.0.4 precision=3.
// Adapted from MIT-licensed PyLightcurve; see ../pylightcurve_LICENSE.txt.
// Retains its tabulated nodes, contact cases, clamps and flux normalization.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>

namespace exotic_claret {
const double pi=3.141592653589793238462643383279502884;
inline double clip(double x,double lo,double hi) { return std::isnan(x)?x:std::max(lo,std::min(hi,x)); }
inline double primitive(const double* a,double r) {
    const double mu4=1.-r*r, mu2=std::sqrt(mu4), mu1=std::sqrt(mu2);
    return -(2.*(1.-a[0]-a[1]-a[2]-a[3])/4.)*mu4
        -(2.*a[0]/5.)*mu4*mu1 -(2.*a[1]/6.)*mu4*mu2
        -(2.*a[2]/7.)*mu4*mu2*mu1 -(2.*a[3]/8.)*mu4*mu4;
}
const double weights[30] = {0.1028526528935588,0.1028526528935588,0.1017623897484055,0.1017623897484055,0.0995934205867953,0.0995934205867953,0.0963687371746443,0.0963687371746443,0.0921225222377861,0.0921225222377861,0.086899787201083,0.086899787201083,0.0807558952294202,0.0807558952294202,0.0737559747377052,0.0737559747377052,0.0659742298821805,0.0659742298821805,0.0574931562176191,0.0574931562176191,0.0484026728305941,0.0484026728305941,0.0387991925696271,0.0387991925696271,0.0287847078833234,0.0287847078833234,0.018466468311091,0.018466468311091,0.0079681924961666,0.0079681924961666};
const double nodes[30] = {-0.0514718425553177,0.0514718425553177,-0.1538699136085835,0.1538699136085835,-0.2546369261678899,0.2546369261678899,-0.3527047255308781,0.3527047255308781,-0.4470337695380892,0.4470337695380892,-0.5366241481420199,0.5366241481420199,-0.6205261829892429,0.6205261829892429,-0.6978504947933158,0.6978504947933158,-0.7677774321048262,0.7677774321048262,-0.8295657623827684,0.8295657623827684,-0.8825605357920527,0.8825605357920527,-0.9262000474292743,0.9262000474292743,-0.9600218649683075,0.9600218649683075,-0.9836681232797472,0.9836681232797472,-0.9968934840746495,0.9968934840746495};

inline double quadrature(const double* a,double radius,double z,double lower,double upper) {
    const double half=(upper-lower)/2., center=(upper+lower)/2.;
    double sum=0.;
    for (int k=0;k<30;++k) {
        const double r=half*nodes[k]+center, rsq=r*r;
        const double mu4=1.-rsq, mu2=std::sqrt(mu4), mu1=std::sqrt(mu2);
        const double intensity=(1.-a[0]-a[1]-a[2]-a[3])+a[0]*mu1+a[1]*mu2+a[2]*mu2*mu1+a[3]*mu4;
        const double argument=(-radius*radius+z*z+rsq)/(2.*z*r);
        const double angle=std::acos(std::isnan(argument)?argument:std::min(argument,1.));
        sum += weights[k]*(intensity*r*angle);
    }
    return half*sum;
}
inline double sector(const double* a,double radius,double z,double theta1,double theta2,bool plus) {
    const double sign=plus?1.:-1.;
    const double s1=z*std::sin(theta1),s2=z*std::sin(theta2);
    const double rr1=clip(z*std::cos(theta1)+sign*std::sqrt(std::max(radius*radius-s1*s1,0.)),0.,1.);
    const double rr2=clip(z*std::cos(theta2)+sign*std::sqrt(std::max(radius*radius-s2*s2,0.)),0.,1.);
    const double w1=std::min(theta1,theta2),w2=std::max(theta1,theta2);
    const double r1=std::min(rr1,rr2),r2=std::max(rr1,rr2);
    const double parta=primitive(a,0.)*(w1-w2);
    const double partb=primitive(a,r1)*(plus?w2:-w1);
    const double partc=primitive(a,r2)*(plus?-w1:w2);
    const double partd=quadrature(a,radius,z,r1,r2);
    return plus?parta+partb+partc+partd:parta+partb+partc-partd;
}
inline double flux(double input_z,double radius,const double* a) {
    double z=input_z<0.?1.+100.*radius:input_z;
    if (!std::isnan(z)) z=std::max(z,1.e-10);
    const double zsq=z*z, sum=z+radius,difference=radius-z, square_difference=zsq-radius*radius;
    const bool c1=z<radius && sum<=1.;
    const bool ca=z<radius && sum>1. && difference<1.;
    const bool cb=z<radius && sum>1. && difference>1.;
    const bool c2=z==radius && sum<=1.;
    const bool cc=z==radius && sum>1.;
    const bool c3=z>radius && sum<1.;
    const bool c4=z>radius && sum==1.;
    const bool c5=z>radius && sum>1. && square_difference<1.;
    const bool c6=z>radius && sum>1. && square_difference==1.;
    const bool c7=z>radius && sum>1. && square_difference>1. && -1.<difference;
    const double ph=std::acos(clip((1.-radius*radius+zsq)/(2.*z),-1.,1.));
    double t1=(c5||ca||cc)?ph:0.;
    const double ratio=radius/z;
    double t2=std::asin(std::isnan(ratio)?ratio:std::min(ratio,1.));
    if (c1||ca) t2=pi;
    if (c2||cc) t2=pi/2.;
    if (c7) t2=ph;
    double upper=0.,lower=0.,star=0.;
    if (c1||c2||c3||c4||c5||ca||cc) upper=sector(a,radius,z,t1,t2,true);
    if (cb) upper=(primitive(a,1.)-primitive(a,0.))*pi;
    if (c3||c4||c5||c6||c7) lower=sector(a,radius,z,0.,t2,false);
    if (c5||c6||c7||ca||cc) star=(primitive(a,1.)-primitive(a,0.))*std::abs(ph);
    const double total=(primitive(a,1.)-primitive(a,0.))*(2.*pi);
    return 1.-(2./total)*(upper+star-lower);
}
} // namespace exotic_claret
