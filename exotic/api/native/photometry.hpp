// EXOTIC native photometry. Double precision; no fast-math.
// Circle overlap translated from Photutils 3.1.1 (Thomas Robitaille), BSD-3.
// Dense TRF and finite differences translated from SciPy 1.18.1, BSD-3.
// See ../photutils_LICENSE.txt and ../scipy_photometry_LICENSE.txt.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
namespace exotic_photo {
using Vec = std::array<double, 7>;
using Matrix = std::vector<double>;
const double nan = std::numeric_limits<double>::quiet_NaN();
const double inf = std::numeric_limits<double>::infinity();
const double eps = std::numeric_limits<double>::epsilon();

// NumPy's contiguous pairwise summation order, including the 128-item block.
inline double sum(const double* p, size_t n) {
    if (n < 8) { double s=-0.0; for(size_t i=0;i<n;++i) s+=p[i]; return s; }
    if (n <= 128) {
        double r[8]; for(int j=0;j<8;++j) r[j]=p[j];
        size_t i=8; for(;i<n-(n%8);i+=8) for(int j=0;j<8;++j) r[j]+=p[i+j];
        double s=((r[0]+r[1])+(r[2]+r[3]))+((r[4]+r[5])+(r[6]+r[7]));
        for(;i<n;++i) s+=p[i]; return s;
    }
    size_t k=n/2; k-=k%8; return sum(p,k)+sum(p+k,n-k);
}
inline double sum(const Matrix& a) { return sum(a.data(),a.size()); }
inline double dot(const Vec& a,const Vec& b) { double s=0;for(int j=0;j<7;++j)s+=a[j]*b[j];return s; }
inline double norm(const Vec& a) { return std::sqrt(dot(a,a)); }
inline Vec mul(const Vec& a,const Vec& b) {Vec c;for(int j=0;j<7;++j)c[j]=a[j]*b[j];return c;}
inline Vec scale(Vec a,double t) {for(double& v:a)v*=t;return a;}
inline Vec add(Vec a,const Vec& b) {for(int j=0;j<7;++j)a[j]+=b[j];return a;}
inline double floor_sqrt(double a) {return a>0?std::sqrt(a):0;}
inline double triangle(double x1,double y1,double x2,double y2,double x3,double y3) {
    return .5*std::abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
}
inline double arc(double x1,double y1,double x2,double y2,double r) {
    double dx=x2-x1,dy=y2-y1,theta=2*std::asin(.5*std::sqrt(dx*dx+dy*dy)/r);
    return .5*r*r*(theta-std::sin(theta));
}
inline double overlap_core(double x0,double y0,double x1,double y1,double r) {
    if(x0*x0+y0*y0>r*r)return 0;
    if(x1*x1+y1*y1<r*r)return (x1-x0)*(y1-y0);
    double d1=floor_sqrt(x1*x1+y0*y0),d2=floor_sqrt(x0*x0+y1*y1),a,b,c,d;
    if(d1<r&&d2<r){a=floor_sqrt(r*r-y1*y1);b=y1;c=x1;d=floor_sqrt(r*r-x1*x1);
        return (x1-x0)*(y1-y0)-triangle(a,b,c,d,x1,y1)+arc(a,b,c,d,r);}
    if(d1<r){a=x0;b=floor_sqrt(r*r-x0*x0);c=x1;d=floor_sqrt(r*r-x1*x1);
        return arc(a,b,c,d,r)+triangle(a,b,a,y0,x1,y0)+triangle(a,b,c,y0,c,d);}
    if(d2<r){a=floor_sqrt(r*r-y0*y0);b=y0;c=floor_sqrt(r*r-y1*y1);d=y1;
        return arc(a,b,c,d,r)+triangle(a,b,x0,b,x0,y1)+triangle(a,b,x0,d,c,d);}
    a=floor_sqrt(r*r-y0*y0);b=y0;c=x0;d=floor_sqrt(r*r-x0*x0);
    return arc(a,b,c,d,r)+triangle(a,b,c,d,x0,y0);
}
inline double overlap(double x0,double y0,double x1,double y1,double r) {
    if(x0>=0){if(y0>=0)return overlap_core(x0,y0,x1,y1,r);
        if(y1<=0)return overlap_core(-y1,x0,-y0,x1,r);
        return overlap(x0,y0,x1,0,r)+overlap(x0,0,x1,y1,r);}
    if(x1<=0){if(y0>=0)return overlap_core(-x1,y0,-x0,y1,r);
        if(y1<=0)return overlap_core(-x1,-y1,-x0,-y0,r);
        return overlap(x0,y0,x1,0,r)+overlap(x0,0,x1,y1,r);}
    if(y0>=0||y1<=0)return overlap(x0,y0,0,y1,r)+overlap(0,y0,x1,y1,r);
    return overlap(x0,y0,0,0,r)+overlap(0,y0,x1,0,r)+overlap(x0,0,0,y1,r)+overlap(0,0,x1,y1,r);
}
struct Mask {
    int x0,y0,nx,ny; Matrix weights;
    Mask(double xc,double yc,double r,double inner,bool center) {
        x0=int(std::floor(xc-r+.5));y0=int(std::floor(yc-r+.5));
        nx=int(std::ceil(xc+r+.5))-x0;ny=int(std::ceil(yc+r+.5))-y0;
        weights.resize(size_t(nx)*ny);
        double xmin=x0-.5-xc,ymin=y0-.5-yc;
        double dx=((x0+nx-.5-xc)-xmin)/nx,dy=((y0+ny-.5-yc)-ymin)/ny;
        auto pixel=[&](int x,int y,double radius){
            double a=xmin+x*dx,b=ymin+y*dy,c=a+dx,d=b+dy;
            double xx=a+dx*.5,yy=b+dy*.5,dist=std::sqrt(xx*xx+yy*yy);
            if(center)return xx*xx+yy*yy<radius*radius?1.:0.;
            double pr=.5*std::sqrt(dx*dx+dy*dy);
            if(dist<radius-pr)return 1.;
            if(dist<radius+pr)return overlap(a,b,c,d,radius)/(dx*dy);
            return 0.;
        };
        for(int y=0;y<ny;++y)for(int x=0;x<nx;++x)
            weights[size_t(y)*nx+x]=pixel(x,y,r)-(inner>0?pixel(x,y,inner):0);
    }
};
struct Pixel {double value,weight;};
inline double percentile(const std::vector<Pixel>& p,double q) {
    if(p.empty())return nan;
    Matrix w;w.reserve(p.size());for(auto& a:p)w.push_back(a.weight);
    double total=sum(w),cum=0,prev=0,prevvalue=p[0].value;
    if(!(total>0)||!std::isfinite(total))return nan;
    for(size_t i=0;i<p.size();++i){cum+=p[i].weight;double c=(cum-.5*p[i].weight)/total;
        if(q<c){if(i==0)return p[0].value;return prevvalue+(q-prev)*(p[i].value-prevvalue)/(c-prev);}
        // np.interp selects the rightmost entry at equal CDF coordinates.
        // Tiny fractional edge weights can round to the same cumulative value.
        prev=c;prevvalue=p[i].value;}
    return p.back().value;
}
inline double scatter(const std::vector<Pixel>& p) {
    Matrix w,v;w.reserve(p.size());v.reserve(p.size());
    for(auto& a:p){w.push_back(a.weight);v.push_back(a.weight*a.value);}
    double total=sum(w),mean=sum(v)/total;
    for(size_t i=0;i<p.size();++i){double d=p[i].value-mean;v[i]=p[i].weight*d*d;}
    return std::sqrt(std::max(sum(v)/total,0.));
}
inline void background(const double* image,int nx,int ny,int ox,int oy,double xc,double yc,
                       double inner,double outer,bool center,double ptol,double sigma,int iters,double* out) {
    Mask mask(xc,yc,outer,inner,center);std::vector<Pixel> p;
    for(int y=0;y<mask.ny;++y)for(int x=0;x<mask.nx;++x){
        int ix=mask.x0+x-ox,iy=mask.y0+y-oy;double w=mask.weights[size_t(y)*mask.nx+x];
        if(ix<0||ix>=nx||iy<0||iy>=ny||!(w>0)||!std::isfinite(w))continue;
        double v=image[size_t(iy)*nx+ix];if(std::isfinite(v))p.push_back({v,w});
    }
    Matrix weights;for(auto& a:p)weights.push_back(a.weight);
    out[2]=sum(weights);out[0]=out[1]=nan;if(p.empty())return;
    if(std::isnan(ptol)){out[2]=0;return;}
    // Keep the original raster order for weighted variance/sums. Sort once for
    // the percentiles; stable equal-valued pixels retain their original order.
    std::vector<Pixel> sorted=p;
    std::stable_sort(sorted.begin(),sorted.end(),[](const Pixel& a,const Pixel& b){return a.value<b.value;});
    double cutoff=percentile(sorted,std::max(0.,std::min(100.,ptol))/100.);
    auto trim=[&](std::vector<Pixel>& v,double limit){v.erase(std::remove_if(v.begin(),v.end(),[&](const Pixel& a){return !(a.value<=limit);}),v.end());};
    trim(p,cutoff);trim(sorted,cutoff);
    for(int i=0;i<iters&&!p.empty();++i){double med=percentile(sorted,.5),sd=scatter(p);
        if(!std::isfinite(med))return;if(!std::isfinite(sd)||sd<=0)break;
        double limit=med+sigma*sd;size_t before=p.size();trim(p,limit);trim(sorted,limit);if(p.size()==before)break;}
    if(!p.empty()){out[0]=percentile(sorted,.5);out[1]=scatter(p);}
}
inline void aperture(const double* image,int nx,int ny,int ox,int oy,double xc,double yc,double r,bool center,double* out) {
    Mask mask(xc,yc,r,0,center);Matrix products(mask.weights.size());
    bool intersects=mask.x0<nx+ox&&mask.x0+mask.nx>ox&&mask.y0<ny+oy&&mask.y0+mask.ny>oy;
    for(int y=0;y<mask.ny;++y)for(int x=0;x<mask.nx;++x){
        int ix=mask.x0+x-ox,iy=mask.y0+y-oy;size_t i=size_t(y)*mask.nx+x;
        products[i]=mask.weights[i]*(ix<0||ix>=nx||iy<0||iy>=ny?0.:image[size_t(iy)*nx+ix]);
    }
    out[0]=intersects?sum(products):0.;out[1]=intersects?sum(mask.weights):nan;out[2]=intersects?1.:0.;
}
// ---- Seven-parameter bounded Gaussian least squares (linear loss, unit scale).
// LAPACK SVD is supplied by the existing SciPy runtime, avoiding a second
// linear-algebra library and preserving the reference decomposition algorithm.
using Svd = int (*)(size_t,const double*,const double*,double*,double*,double*);
inline Matrix residual(const double* xs,const double* ys,const double* data,size_t m,const Vec& p) {
    Matrix f(m);double c=std::cos(p[5]),s=std::sin(p[5]);
    for(size_t i=0;i<m;++i){double dx=xs[i]-p[0],dy=ys[i]-p[1],rx=dx*c-dy*s,ry=dx*s+dy*c;
        double gx=std::exp(-rx*rx/(2*p[3]*p[3])),gy=std::exp(-ry*ry/(2*p[4]*p[4]));
        f[i]=data[i]-(p[2]*gx*gy+p[6]);}return f;
}
inline Matrix jacobian(const double* xs,const double* ys,const double* data,size_t m,const Vec& p,
                       const Matrix& f,const Vec& lo,const Vec& hi,int method) {
    Matrix J(m*7);double rel=std::pow(eps,method==2?.5:1./3.);
    for(int j=0;j<7;++j){double h=rel*(p[j]>=0?1.:-1.)*std::max(1.,std::abs(p[j]));
        double ld=p[j]-lo[j],ud=hi[j]-p[j];bool one=method==2;
        if(method==2){bool violated=p[j]+h<lo[j]||p[j]+h>hi[j],fitting=std::abs(h)<=std::max(ld,ud);
            if(violated&&fitting)h=-h;if(!fitting)h=ud>=ld?ud:-ld;
        }else{h=std::abs(h);if(!(ld>=h&&ud>=h)){one=true;h=ud>=ld?std::min(h,.5*ud):-std::min(h,.5*ld);
                double md=std::min(ud,ld);if(std::abs(h)<=md){h=md;one=false;}}}
        Vec a=p,b=p;a[j]=p[j]+(method==3&&!one?-h:h);b[j]=p[j]+(one?2*h:h);
        Matrix fa=residual(xs,ys,data,m,a),fb;if(method==3)fb=residual(xs,ys,data,m,b);
        double dx=method==2?a[j]-p[j]:(one?b[j]-p[j]:b[j]-a[j]);
        for(size_t i=0;i<m;++i)J[i*7+j]=(method==2?fa[i]-f[i]:(one?-3.*f[i]+4.*fa[i]-fb[i]:fb[i]-fa[i]))/dx;
    }return J;
}
inline Vec gradient(const Matrix& J,const Matrix& f) {
    Vec g{};for(size_t i=0;i<f.size();++i)for(int j=0;j<7;++j)g[j]+=J[i*7+j]*f[i];return g;
}
inline void cl_scale(const Vec& x,const Vec& g,const Vec& lo,const Vec& hi,Vec& v,Vec& dv) {
    for(int j=0;j<7;++j){v[j]=1;dv[j]=0;if(g[j]<0&&std::isfinite(hi[j])){v[j]=hi[j]-x[j];dv[j]=-1;}
        if(g[j]>0&&std::isfinite(lo[j])){v[j]=x[j]-lo[j];dv[j]=1;}}
}
inline Vec feasible(Vec x,const Vec& lo,const Vec& hi,double step) {
    for(int j=0;j<7;++j){double ld=x[j]-lo[j],ud=hi[j]-x[j];int active=0;
        if(std::isfinite(lo[j])&&ld<=std::min(ud,step*std::max(1.,std::abs(lo[j]))))active=-1;
        if(std::isfinite(hi[j])&&ud<=std::min(ld,step*std::max(1.,std::abs(hi[j]))))active=1;
        if(active<0)x[j]=step?lo[j]+step*std::max(1.,std::abs(lo[j])):std::nextafter(lo[j],hi[j]);
        if(active>0)x[j]=step?hi[j]-step*std::max(1.,std::abs(hi[j])):std::nextafter(hi[j],lo[j]);
        if(x[j]<lo[j]||x[j]>hi[j])x[j]=.5*(lo[j]+hi[j]);}return x;
}
inline double bound_stride(const Vec& x,const Vec& s,const Vec& lo,const Vec& hi,Vec* hits=nullptr) {
    Vec steps;double t=inf;for(int j=0;j<7;++j){steps[j]=s[j]!=0?std::max((lo[j]-x[j])/s[j],(hi[j]-x[j])/s[j]):inf;t=std::min(t,steps[j]);}
    if(hits)for(int j=0;j<7;++j)(*hits)[j]=steps[j]==t?(s[j]>0?1:-1):0;return t;
}
inline void quadratic(const Matrix& J,const Vec& g,const Vec& s,const Vec& diag,double& a,double& b,double& c,const Vec* origin=nullptr) {
    double vv=0,uv=0,uu=0;for(size_t i=0;i<J.size()/7;++i){double v=0,u=0;for(int j=0;j<7;++j){v+=J[i*7+j]*s[j];if(origin)u+=J[i*7+j]*(*origin)[j];}vv+=v*v;uv+=u*v;uu+=u*u;}
    a=.5*(vv+dot(mul(s,diag),s));b=dot(g,s);c=0;
    if(origin){b+=uv+dot(mul(*origin,diag),s);c=.5*uu+dot(g,*origin)+.5*dot(mul(*origin,diag),*origin);}
}
inline double quad_value(const Matrix& J,const Vec& g,const Vec& s,const Vec& diag) {double a,b,c;quadratic(J,g,s,diag,a,b,c);return a+b;}
inline double minimize(double a,double b,double l,double u,double c,double& value) {
    double t=l;value=l*(a*l+b)+c;double v=u*(a*u+b)+c;if(v<value){value=v;t=u;}
    if(a!=0){double z=-.5*b/a;if(l<z&&z<u){v=z*(a*z+b)+c;if(v<value){value=v;t=z;}}}return t;
}
inline Vec select(const Vec& x,const Matrix& J,const Vec& diag,const Vec& g,Vec p,Vec ph,const Vec& d,
                  double delta,const Vec& lo,const Vec& hi,double theta,Vec& sh,double& predicted) {
    Vec xp=add(x,p);bool inside=true;for(int j=0;j<7;++j)inside=inside&&xp[j]>=lo[j]&&xp[j]<=hi[j];
    if(inside){sh=ph;predicted=-quad_value(J,g,ph,diag);return p;}
    Vec hits;double stride=bound_stride(x,p,lo,hi,&hits);Vec rh=ph;
    for(int j=0;j<7;++j)if(hits[j])rh[j]*=-1;
    Vec r=mul(d,rh);p=scale(p,stride);ph=scale(ph,stride);Vec on=add(x,p);
    double aa=dot(rh,rh),bb=dot(ph,rh),cc=dot(ph,ph)-delta*delta;
    if(aa==0||cc>0)throw std::runtime_error("Invalid PSF trust-region intersection");
    double q=-(bb+std::copysign(std::sqrt(bb*bb-aa*cc),bb));
    double to_tr=std::max(q/aa,cc/q),to_bound=bound_stride(on,r,lo,hi);
    double rs=std::min(to_bound,to_tr),rl=rs>0?(1-theta)*stride/rs:0,ru=rs>0?(rs==to_bound?theta*to_bound:to_tr):-1;
    double rv=inf,a,b,c;
    if(rl<=ru){quadratic(J,g,rh,diag,a,b,c,&ph);rs=minimize(a,b,rl,ru,c,rv);rh=add(scale(rh,rs),ph);r=mul(rh,d);}
    p=scale(p,theta);ph=scale(ph,theta);double pv=quad_value(J,g,ph,diag);
    Vec ah=scale(g,-1),ag=mul(d,ah);to_tr=delta/norm(ah);to_bound=bound_stride(x,ag,lo,hi);
    double as=to_bound<to_tr?theta*to_bound:to_tr,av;quadratic(J,g,ah,diag,a,b,c);as=minimize(a,b,0,as,0,av);
    ah=scale(ah,as);ag=scale(ag,as);
    if(pv<rv&&pv<av){sh=ph;predicted=-pv;return p;}
    if(rv<pv&&rv<av){sh=rh;predicted=-rv;return r;}
    sh=ah;predicted=-av;return ag;
}
inline Vec trust(const Vec& sv,const Vec& uf,const double* V,double delta,size_t m,double& alpha) {
    Vec suf=mul(sv,uf),p{};bool full=m>=7&&sv[6]>eps*m*sv[0];
    auto transform=[&](const Vec& v){Vec r{};for(int j=0;j<7;++j)for(int k=0;k<7;++k)r[j]-=V[j*7+k]*v[k];return r;};
    if(full){Vec z;for(int j=0;j<7;++j)z[j]=uf[j]/sv[j];p=transform(z);if(norm(p)<=delta){alpha=0;return p;}}
    auto phi=[&](double al,double& derivative){Vec z;double dp=0;for(int j=0;j<7;++j){double den=sv[j]*sv[j]+al;z[j]=suf[j]/den;dp+=suf[j]*suf[j]/(den*den*den);}double pn=norm(z);derivative=-dp/pn;return pn-delta;};
    double upper=norm(suf)/delta,lower=0,derivative;
    if(full){double f=phi(0,derivative);lower=-f/derivative;}
    if(!full&&alpha==0)alpha=std::max(.001*upper,std::sqrt(lower*upper));
    for(int i=0;i<10;++i){if(alpha<lower||alpha>upper)alpha=std::max(.001*upper,std::sqrt(lower*upper));
        double f=phi(alpha,derivative);if(f<0)upper=alpha;double ratio=f/derivative;lower=std::max(lower,alpha-ratio);alpha-=(f+delta)*ratio/delta;if(std::abs(f)<.01*delta)break;}
    Vec z;for(int j=0;j<7;++j)z[j]=suf[j]/(sv[j]*sv[j]+alpha);p=transform(z);return scale(p,delta/norm(p));
}
inline void fit(const double* xs,const double* ys,const double* data,size_t m,const double* initial,
                const double* lower,const double* upper,int method,Svd svd,double* output,double* info) {
    Vec x,lo,hi;for(int j=0;j<7;++j){x[j]=initial[j];lo[j]=lower[j];hi[j]=upper[j];
        if(!std::isfinite(x[j])||!(lo[j]<hi[j])||x[j]<lo[j]||x[j]>hi[j])throw std::invalid_argument("Invalid PSF initial point or bounds");}
    x=feasible(x,lo,hi,1e-10);Matrix f=residual(xs,ys,data,m,x);for(double v:f)if(!std::isfinite(v))throw std::invalid_argument("PSF residuals are not finite at the initial point");
    Matrix J=jacobian(xs,ys,data,m,x,f,lo,hi,method);Vec g=gradient(J,f),v,dv;cl_scale(x,g,lo,hi,v,dv);
    Vec scaled;for(int j=0;j<7;++j)scaled[j]=x[j]/std::sqrt(v[j]);double delta=norm(scaled);if(delta==0)delta=1;
    auto cost_of=[](const Matrix& ff){double c=0;for(double a:ff)c+=a*a;return .5*c;};
    double cost=cost_of(f),alpha=0;int nfev=1,njev=1,status=0;double gnorm=inf;
    // Reference SciPy default max_nfev = 100 * number of parameters.
    while(true){cl_scale(x,g,lo,hi,v,dv);gnorm=0;for(int j=0;j<7;++j)gnorm=std::max(gnorm,std::abs(g[j]*v[j]));
        if(gnorm<1e-8)status=1;if(status||nfev==700)break;
        Vec d,diag,gh;for(int j=0;j<7;++j){d[j]=std::sqrt(v[j]);diag[j]=g[j]*dv[j];gh[j]=d[j]*g[j];}
        Matrix Jh(m*7),aug((m+7)*7,0),fa(m+7,0);for(size_t i=0;i<m;++i){fa[i]=f[i];for(int j=0;j<7;++j)aug[i*7+j]=Jh[i*7+j]=J[i*7+j]*d[j];}
        for(int j=0;j<7;++j)aug[(m+j)*7+j]=std::sqrt(diag[j]);
        Vec sv,uf;double V[49];if(svd(m+7,aug.data(),fa.data(),sv.data(),V,uf.data()))throw std::runtime_error("PSF SVD failed");
        double theta=std::max(.995,1-gnorm),reduction=-1;Vec xn;Matrix fn;double costnew=cost;
        while(reduction<=0&&nfev<700){Vec ph=trust(sv,uf,V,delta,m,alpha),p=mul(d,ph),sh;double predicted;
            Vec step=select(x,Jh,diag,gh,p,ph,d,delta,lo,hi,theta,sh,predicted);xn=feasible(add(x,step),lo,hi,0);
            fn=residual(xs,ys,data,m,xn);++nfev;double shnorm=norm(sh);bool finite=true;for(double a:fn)finite=finite&&std::isfinite(a);
            if(!finite){delta=.25*shnorm;continue;}
            costnew=cost_of(fn);reduction=cost-costnew;double ratio=predicted>0?reduction/predicted:(predicted==0&&reduction==0?1:0),next=delta;
            if(ratio<.25)next=.25*shnorm;else if(ratio>.75&&shnorm>.95*delta)next*=2;
            bool ft=reduction<1e-8*cost&&ratio>.25,xt=norm(step)<1e-8*(1e-8+norm(x));status=ft?(xt?4:2):(xt?3:0);
            if(status)break;alpha*=delta/next;delta=next;
        }
        if(reduction>0){x=xn;f=fn;cost=costnew;J=jacobian(xs,ys,data,m,x,f,lo,hi,method);++njev;g=gradient(J,f);}
    }
    for(int j=0;j<7;++j)output[j]=x[j];info[0]=cost;info[1]=nfev;info[2]=njev;info[3]=status;info[4]=gnorm;
}
} // namespace exotic_photo
