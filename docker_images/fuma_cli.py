from python_on_whales import docker


def run_find_rsid(
	container_name: str,
	ref_data_path_on_host: str,
	jobs_dir_on_host: str,
	image_name: str,
	job_location: str,
) -> None:
	docker.run(
		image_name,
		[
			"/bin/sh",
			"-c",
			(
				f"python find_rsid_from_grch38.py --filedir {job_location}/ "
				f">>{job_location}/job.log 2>>{job_location}/error.log"
			),
		],
		remove=True,
		networks=["none"],
		name=container_name,
		volumes=[
			(ref_data_path_on_host, "/data"),
			(jobs_dir_on_host, jobs_dir_on_host),
		],
	)
 
def gwas_file(
	container_name: str,
	ref_data_path_on_host: str,
	jobs_dir_on_host: str,
	image_name: str,
	job_location: str,
) -> None:
	docker.run(
		image_name,
		[
			"/bin/sh",
			"-c",
			(
				f"python gwas_file.py {job_location}/ "
				f">>{job_location}/job.log 2>>{job_location}/error.log"
			),
		],
		remove=True,
		networks=["none"],
		name=container_name,
		volumes=[
			(ref_data_path_on_host, "/data"),
			(jobs_dir_on_host, jobs_dir_on_host),
		],
	)
 
# run_find_rsid(
#     container_name="find_rsid_container",
#     ref_data_path_on_host="/mnt/d/reference_data",
#     jobs_dir_on_host="/fuma_data/users_data",
#     image_name="docker_images-gwas_file_cli ",
#     job_location="/fuma_data/users_data/jobs/1"
# )

gwas_file(
    container_name="gwas_file_container",
    ref_data_path_on_host="/mnt/d/reference_data",
    jobs_dir_on_host="/fuma_data/users_data",
    image_name="docker_images-gwas_file_cli",
    job_location="/fuma_data/users_data/jobs/1")