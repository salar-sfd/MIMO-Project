function detector = makeDetector(method, cons, consEnergy, modulation)
    [methodName, params] = parseMethodString(method);

    switch upper(methodName)
        case 'ZF'
            detector = @(y_v,H_m,snr,N0,Nt,Nr) ZF(y_v, H_m, cons);

        case 'MMSE'
            detector = @(y_v,H_m,snr,N0,Nt,Nr) MMSE(y_v, H_m, snr, cons);

        case 'LRA'
            alpha = getOr(params, 'alpha', 0.75);
            method = getOr(params, 'method', 'ZF');
            detector = @(y_v, H_m, snr, N0, Nt, Nr) LRA(y_v, H_m, snr, cons, consEnergy, alpha, method, modulation);
            
        case 'SD'
            d = getOr(params,'d',0.7);
            detector = @(y_v,H_m,snr,N0,Nt,Nr) SD(y_v, H_m, d, cons, consEnergy, modulation);

        case 'OGD'
            detector = @(y_v,H_m,snr,N0,Nt,Nr) OGD(y_v, H_m, cons, consEnergy, modulation);

        otherwise
            error('Unknown detection method: %s', methodName);
    end
end

function v = getOr(params, field, default)
    if isfield(params, field)
        v = params.(field);
    else
        v = default;
    end
end

function [name, params] = parseMethodString(s)
    params = struct();
    if isstring(s), s = char(s); end
    s = strtrim(s);
    if numel(s)>=2 && ((s(1)=='"' && s(end)=='"')||(s(1)=='''' && s(end)=='''')), s=s(2:end-1); s=strtrim(s); end

    name = strtrim(regexp(s, '^[^(]+', 'match', 'once'));
    argCell = regexp(s, '\((.*?)\)', 'tokens', 'once'); 
    if isempty(argCell), return; end
    argStr = argCell{1};
    parts = strsplit(argStr, ',');
    for i=1:numel(parts)
        kv = strtrim(parts{i});
        if isempty(kv), continue; end
        idx = regexp(kv, '=', 'once');
        if isempty(idx)
            key = matlab.lang.makeValidName(kv);
            params.(key) = true;
        else
            key = matlab.lang.makeValidName(strtrim(kv(1:idx-1)));
            val = strtrim(kv(idx+1:end));
            if numel(val)>=2 && ((val(1)=='"'&&val(end)=='"')||(val(1)==''''&&val(end)=='''')), val=val(2:end-1); params.(key)=val; continue; end
            vnum = str2double(val);
            if ~isnan(vnum), params.(key)=vnum; else params.(key)=val; end
        end
    end
end

