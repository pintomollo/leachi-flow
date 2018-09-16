## Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
## 
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
## 
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
## 
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} function_name ()
## @end deftypefn

function varargout = subsref (obj, idx)

  persistent __method__ method4field typeNotImplemented
  if isempty(__method__)

    __method__ = struct();

    __method__.getPosition = @(o,varargin) getPosition (o, varargin{:});
    __method__.setPosition = @(o,varargin) setPosition (o, varargin{:});
    __method__.setClosed = @(o,varargin) setClosed (o, varargin{:});
    __method__.delete = @(o,varargin) delete (o, varargin{:});
    __method__.display = @(o,varargin) display (o, varargin{:});

    # Error strings
    method4field = "Class #s has no field #s. Use #s() for the method.";
    typeNotImplemented = "#s no implemented for class #s.";

  end

  if ( !strcmp (class (obj), 'impoly') )
    error ("Object must be of the impoly class but '#s' was used", class (obj) );
  elseif ( idx(1).type != '.' )
    error ("Invalid index for class #s", class (obj) );
  endif

  method = idx(1).subs;
  if ~isfield(__method__, method)
    error('Unknown method #s.',method);
  else
    fhandle = __method__.(method);
  end

  if numel (idx) == 1 # can't access properties, only methods

    error (method4field, class (obj), method, method);

  end

  if strcmp (idx(2).type, '()')

    args = idx(2).subs;

    if (method(1) == 'g')
      if isempty(args)
        out = fhandle (obj);
      else
        out = fhandle (obj, args{:});
      end
      varargout{1} = out;
    else

      if isempty(args)
        fhandle (obj);
      else
        fhandle (obj, args{:});
      end
    end

  else

    error (typeNotImplemented,[method idx(2).type], class (obj));

  end

endfunction
