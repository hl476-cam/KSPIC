function [wfun,tune] = statrobustwfun_IC(wfun,tune,alpha)

% Modified based on STATROBUSTWFUN.

%STATROBUSTWFUN Get robust weighting function and tuning constant

% Copyright 2005-2007 The MathWorks, Inc.


% Convert name of weight function to a handle to a local function, and get
% the default value of the tuning parameter
t = [];
if ischar(wfun)
    switch(wfun)
        case 'andrews'
            wfun = @andrews;
            t = 1.339;
        case 'bisquare'
            wfun = @bisquare;
            t = 4.685;
        case 'cauchy'
            wfun = @cauchy;
            t= 2.385;
        case 'fair'
            wfun = @fair;
            t = 1.400;
        case 'huber'
            wfun = @huber;
            t = 1.345;
        case 'logistic'
            wfun = @logistic;
            t = 1.205;
        case 'ols'
            wfun = @ols;
            t = 1;
        case 'talwar'
            wfun = @talwar;
            t = 2.795;
        case 'welsch'
            wfun = @welsch;
            t = 2.985;
        case 'MRA'
            wfun = @MRA;
            t = 2.985;
        case 'MRA2'
            wfun = @MRA2;
            t = 2.985;
        case 'MRA3'
            wfun = @MRA3;
            t = 2.985;
        case 'MRA4'
            wfun = @MRA4;
            t = 2.985;
        case 'MRA5'
            wfun = @MRA5;
            t = 2.985;

    end
end

% Use the default tuning parameter or check the supplied one
if isempty(tune)
    if isempty(t)
        tune = 1;
    else
        tune = t;
    end
elseif (tune<=0)
    m = message('stats:statrobustwfun:BadTuningConstant');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end

% --------- weight functions

function w = andrews(r)
r = max(sqrt(eps(class(r))), abs(r));
w = (abs(r)<pi) .* sin(r) ./ r;

function w = bisquare(r)
w = (abs(r)<1) .* (1 - r.^2).^2;

function w = cauchy(r)
w = 1 ./ (1 + r.^2);

function w = fair(r)
w = 1 ./ (1 + abs(r));

function w = huber(r)
w = 1 ./ max(1, abs(r));

function w = logistic(r)
r = max(sqrt(eps(class(r))), abs(r));
w = tanh(r) ./ r;

function w = ols(r)
w = ones(size(r));

function w = talwar(r)
w = 1 * (abs(r)<1);

function w = welsch(r)
w = exp(-(r.^2));

function w = MRA(r,R)
% w = exp(-((r/1).^2))'*(20-20*exp(-((R/100).^2)));
% w = exp(-((r/1).^2))'*R;
w = exp(-(r.^2))'.*exp(-((1-R).^2));
% w = exp(-(r.^2))'*R/R;


function w = MRA2(r,X)
w = exp(-(r.^2)).*abs(X);


function w = MRA3(theta,R)
w = exp(-(theta.^2)).*abs(R);

function w = MRA4(theta)
w = exp(-(theta.^2));

function w = MRA5(theta,x,y,alpha)
% w = exp(-(theta.^2));
wa = exp(-(theta.^2));
wv=((abs(x)+abs(y)/2).^2).^alpha;
w=wv.*wa;

