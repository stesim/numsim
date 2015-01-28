function visualize(ranks,step,instance,varargin)
	funnames = {'u','v','p','uv','T','zeta'};

	if length(varargin) == 0
		varargin = funnames;
	end

	isset = @(funname)(any(strcmp(funname,varargin)));

	step_str = num2str(step,'%06d');
	rx = ranks(2);
	ry = ranks(1);

	uvp = stitch_subdomains('uvpT',6);
    
	if isset('u')
		figure(1);
		surf(uvp{1},uvp{2},uvp{3},'EdgeColor','None');
		title('velocity_x');
		xlabel('x');
		ylabel('y');
	end

	if isset('v')
		figure(2);
		surf(uvp{1},uvp{2},uvp{4},'EdgeColor','None');
		title('velocity_y');
		xlabel('x');
		ylabel('y');
	end

	if isset('p')
		figure(3);
		surf(uvp{1},uvp{2},uvp{5},'EdgeColor','None');
		title('pressure');
		xlabel('x');
		ylabel('y');
	end
	
	if isset('uv')
		figure(4);
		clf(4);
		startx = rand([1,floor(gridSize(2)*1.7)]);
		streamslice(uvp{1},uvp{2},uvp{3},uvp{4},1.42);
		title('streamlines');
		xlabel('x');
		ylabel('y');
	end

	%psi = stitch_subdomains('psi',3);

	%figure(5);
	%surf(psi{1},psi{2},psi{3},'EdgeColor','None');
	%title('velocity potential');
	%xlabel('x');
	%ylabel('y');
    
	if isset('T')
		figure(5);
		surf(uvp{1},uvp{2},uvp{6},'EdgeColor','None');
		view(2);
		title('temperature');
		xlabel('x');
		ylabel('y');
	end

	if isset('zeta')
		zeta = stitch_subdomains('zeta',3);

		figure(6);
		surf(zeta{1},zeta{2},zeta{3},'EdgeColor','None');
		title('vorticity');
		xlabel('x');
		ylabel('y');
	end

	function fields = stitch_subdomains(field_names,num_fields)
		size_global = [0,0];

		fields = cell(1,num_fields);

		scan_mask = '';
		for n=1:num_fields
			scan_mask = [scan_mask '%f '];
		end

		for i=1:rx
			for j=1:ry
				filename = ['field-' field_names '-instance-' num2str(instance,'%03d') '-step-' step_str '-rank-' num2str(i-1,'%03d') '-' num2str(j-1,'%03d') '.txt'];
				file = fopen(filename);

				gridSize = fscanf(file, '%d %d');
				gridSize = gridSize([2,1]);
				if i==1 && j==1
					size_global(1) = ry*gridSize(1);
					size_global(2) = rx*gridSize(2);

					for n=1:num_fields
						fields{n} = zeros(size_global);
					end
				end

				fseek(file,0,-1);
				fgetl(file);

				entries = fscanf(file, scan_mask, [num_fields,Inf]);

				sel_y = (j-1)*gridSize(1)+1:j*gridSize(1);
				sel_x = (i-1)*gridSize(2)+1:i*gridSize(2);

				for n=1:num_fields
					fields{n}(sel_y,sel_x) = reshape(entries(n,:),[gridSize(2),gridSize(1)]).';
				end

				fclose(file);
			end
		end
	end
end
