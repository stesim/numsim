function visualize(ranks,step)
	step_str = num2str(step,'%06d');
	files = dir(['field-uvp-step-' step_str '-rank-*-*.txt']);

	size_global = [0;0];
	x_global = zeros(0);
	y_global = zeros(0);
	u_global = zeros(0);
	v_global = zeros(0);
	p_global = zeros(0);

	rx = ranks(1);
	ry = ranks(2);
	for i=1:rx
		for j=1:ry
			filename = ['field-uvp-step-' step_str '-rank-' num2str(i-1,'%03d') '-' num2str(j-1,'%03d') '.txt'];
			file = fopen(filename);

			gridSize = fscanf(file, '%d %d');
			gridSize = gridSize(1:2);
			if i==1 && j==1
				size_global(1) = rx*gridSize(1);
				size_global(2) = ry*gridSize(2);

				x_global = zeros(size_global(2),size_global(1));
				y_global = zeros(size_global(2),size_global(1));
				u_global = zeros(size_global(2),size_global(1));
				v_global = zeros(size_global(2),size_global(1));
				p_global = zeros(size_global(2),size_global(1));
			end

			fseek(file,0,-1);
			fgetl(file);
			entries = fscanf(file, '%f %f %f %f %f', [5,Inf]);
			x = reshape(entries(1,:),[gridSize(1),gridSize(2)]).';
			y = reshape(entries(2,:),[gridSize(1),gridSize(2)]).';
			u = reshape(entries(3,:),[gridSize(1),gridSize(2)]).';
			v = reshape(entries(4,:),[gridSize(1),gridSize(2)]).';
			p = reshape(entries(5,:),[gridSize(1),gridSize(2)]).';

			sel_y = (j-1)*gridSize(2)+1:j*gridSize(2);
			sel_x = (i-1)*gridSize(1)+1:i*gridSize(1);
			x_global(sel_y,sel_x) = x;
			y_global(sel_y,sel_x) = y;
			u_global(sel_y,sel_x) = u;
			v_global(sel_y,sel_x) = v;
			p_global(sel_y,sel_x) = p;
		end
	end

	figure(1);
	surf(x_global,y_global,u_global,'EdgeColor','None');
	title('u');
	xlabel('x');
	xlabel('y');

	figure(2);
	surf(x_global,y_global,v_global,'EdgeColor','None');
	title('v');
	xlabel('x');
	xlabel('y');

	figure(3);
	surf(x_global,y_global,p_global,'EdgeColor','None');
	title('p');
	xlabel('x');
	xlabel('y');
	
	figure(4);
	startx = rand([1,floor(gridSize(1)*1.7)]);
	streamline(x_global,y_global,u_global,v_global,startx,rand([1,size(startx)]));
	title('streamlines');

end
