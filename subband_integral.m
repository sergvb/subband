clc, clear all;

K = 400;
N = 16 ;
fd = 16E3 ;
fs = 1.9313E3 ;
w1 = -2*pi*3E3/fd;
w2 = 2*pi*3E3/fd;

s = cos(2*pi*fs/fd*(0:N-1));

res = zeros(1, K);

for i = 1:K
    h = (w2-w1)/i;
    w = w1:h:w2;
    a = 0;
    for k = 1:i
        t = sum(s.*exp(1i*w(k)*(0:N-1)));
        a = a + t*conj(t);
    end
    res(i) = a*h;
end

semilogx(res);
grid on;