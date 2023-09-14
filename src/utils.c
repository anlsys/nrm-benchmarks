/*******************************************************************************
 * Copyright 2021 UChicago Argonne, LLC.
 * (c.f. AUTHORS, LICENSE)
 *
 * This file is part of the NRM Benchmarks project.
 * For more info, see https://github.com/anlsys/nrm-benchmarks
 *
 * SPDX-License-Identifier: BSD-3-Clause
 ******************************************************************************/

#include "config.h"

#include "nrm-benchmarks.h"

int nrmb_check_double(double ref, double value, int bits)
{
	double diff = NRMB_ABS(ref - value);
	return diff <= NRMB_MAX(NRMB_ABS(ref), NRMB_ABS(value)) * DBL_EPSILON *
		((1 << bits) - 1);
}

int nrmb_check_double_prec(double ref, double value, double eps)
{
	double diff = NRMB_ABS(ref - value);
	return diff <= NRMB_MAX(NRMB_ABS(ref), NRMB_ABS(value)) * eps;
}

static nrm_client_t *nrmb_client;
static nrm_sensor_t *nrmb_sensor;
static nrm_scope_t *nrmb_scope;

int nrmb_init(const char *progname)
{

	nrm_init(NULL, NULL);
	nrm_log_init(stderr, argv[0]);
	nrm_client_create(&nrmb_client, nrm_upstream_uri, nrm_upstream_pub_port,
			  nrm_upstream_rpc_port);
	assert(nrmb_client != NULL);
	nrm_string_t progress_name = nrm_string_fromprintf("nrm.benchmarks.%s",
							   argv[0]);
	nrmb_sensor = nrm_sensor_create(progress_name);
	nrm_client_add_sensor(nrmb_client, nrmb_sensor);

	nrm_vector_t *nrmd_scopes;
	size_t numscopes = 0;
	nrm_client_list_scopes(nrmb_client, &nrmd_scopes);
	nrm_vector_length(nrmd_scopes, &numscopes);
	for (size_t i = 0; i < numscopes; i++) {
		nrm_scope_t *s;
		nrm_vector_pop_back(nrmd_scopes, &s);
		if (!strcmp(s->uuid, "nrm.scope.hwloc.Machine.0")) {
			nrmb_scope = s;
			continue;
		}
		nrm_scope_destroy(s);
	}
	nrm_vector_destroy(&nrmd_scopes);
}

int nrmb_finalize(void)
{
	nrm_client_remove_sensor(nrmb_client, nrmb_sensor);
	nrm_sensor_destroy(&nrmb_sensor);
	nrm_scope_destroy(nrmb_scope);
	nrm_client_destroy(&nrmb_client);
}

int nrmb_send_progress(double value)
{
	nrm_time_t now;
	nrm_time_gettime(&now);
	nrm_send_event(nrmb_client, now, nrmb_sensor, nrmb_scope, value);
}
