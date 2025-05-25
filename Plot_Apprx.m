% ________Eropean Call option under Black-Schole-Merton Model
d = 1.0/(sqrt(2*pi));
r=0.03;
tm = 0.2;           % Time to Maturity
p=1;                % Moneyness
pd=p*exp(-r*tm);    % Discounted moneyness
B = 1;              % Domain limit
n1=10; dx=B/n1; M=zeros(n1,n1);
n2=100; z=zeros(n2,n2);sig1=zeros(1,n2);sig2=zeros(1,n2);
for ii=1:n2
    sig1(ii)=0.04+ (ii)/n2;
    for jj=1:n2
        sig2(jj)=0.04+ (jj)/n2;
        for j= 1:n1
           sj = j*dx;
           eta_j=BS(sj,p*sj,tm,sig1(ii),r);
           LHS=(eta_j/sj - (1-pd)/2)/sig1(ii);
                for i=1:n1
                    si=i*dx;
                    eta_i=BS(si,p*si,tm,sig2(jj),r);
                     RHS=(eta_i/si - (1-pd)/2)/sig2(jj);  
                     M(i,j)=abs(LHS-RHS)/max(abs(LHS),abs(RHS));
                end 
        end 
        z(ii,jj)=max(max(M));
    end
end
max(max(z))
    surf (sig2,sig1,z,"EdgeColor","none")
    map = [(1:-0.02:0)'*[1 1],(1:-0.01:0.50)'];
    colormap(map); colorbar
    xlabel('\sigma_2'), ylabel('\sigma_1'), title('Approximation Error')
    view(90,-90);
    xlim([min(sig2) 1]); ylim([min(sig1) 1]);
    grid on
sigDiff=[]; z_linear=[];
for ii =1:n2
    for jj=1:n2
        sigDiff(end+1)=max(ii,jj)/min(ii,jj);
        z_linear(end+1)= z(ii,jj);
    end
end
sc=scatter(sigDiff,z_linear,5,[0.25 0.25 0.65],'o');
h=gca; set(h,'xscale','log'); grid on
xlabel('max(\sigma_1, \sigma_2)/min(\sigma_1, \sigma_2)')
ylabel('Relative Error')

function [eta]=BS(s,st,tm,sig,r)
        rp= r+ 0.5 *sig^2;
        rm= r- 0.5 *sig^2;
        dn= sig*sqrt(tm);% tm is time to expiry:=T-t
        sk = log(s/st);
        term1= normcdf((sk+rp*tm)/dn);
        term2= normcdf((sk+rm*tm)/dn);
        eta= s*term1- st*exp(-r*tm)*term2;
end