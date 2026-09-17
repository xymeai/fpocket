/*
 * Thin C wrapper around fpocket's internal API, providing a simplified
 * interface suitable for FFI bindings (CFFI).
 *
 * Exposes pocket detection via a flat struct array rather than linked lists.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "fpocket.h"
#include "fparams.h"
#include "rpdb.h"
#include "pocket.h"
#include "descriptors.h"
#include "memhandler.h"
#include "read_mmcif.h"

/* Flat pocket result for FFI consumption */
typedef struct {
    int rank;
    float score;
    float druggability_score;
    int num_alpha_spheres;
    int num_apolar_alpha_spheres;
    int num_polar_alpha_spheres;
    float volume;
    float convex_hull_volume;
    float hydrophobicity_score;
    float polarity_score;
    float charge_score;
    float prop_polar_atoms;
    float mean_alpha_sphere_radius;
    float alpha_sphere_density;
    float mean_local_hydrophobic_density;
    float apolar_alpha_sphere_proportion;
    float max_alpha_sphere_distance;
    float flexibility;
    float surf_vdw;
    float surf_pol_vdw;
    float surf_apol_vdw;
    float bary_x;
    float bary_y;
    float bary_z;
    int inter_chain;
    char chain1[256];
    char chain2[256];
    char lig_tag[8];
    int aa_composition[20];
    int n_abpa;
} fpocket_pocket_result;

/* Result container */
typedef struct {
    fpocket_pocket_result *pockets;
    int num_pockets;
    int error_code;
    char error_message[512];
} fpocket_result;


static s_pdb *open_and_read_pdb(const char *pdb_path, int keep_lig,
                                 s_fparams *params) {
    s_pdb *pdb;
    if (strstr(pdb_path, ".cif")) {
        pdb = open_mmcif((char *)pdb_path, NULL, keep_lig,
                         params->model_number, params);
        if (pdb)
            read_mmcif(pdb, NULL, keep_lig, params->model_number, params);
    } else {
        pdb = rpdb_open((char *)pdb_path, NULL, keep_lig,
                        params->model_number, params);
        if (pdb)
            rpdb_read(pdb, NULL, keep_lig, params->model_number, params);
    }
    return pdb;
}


fpocket_result *fpocket_run(const char *pdb_path,
                             float min_alpha_size,
                             float max_alpha_size,
                             int min_spheres_per_pocket,
                             float clustering_distance,
                             int mc_iterations) {
    fpocket_result *result = calloc(1, sizeof(fpocket_result));
    if (!result) return NULL;

    /* Validate input */
    if (!pdb_path || strlen(pdb_path) == 0) {
        result->error_code = 1;
        snprintf(result->error_message, sizeof(result->error_message),
                 "PDB path is empty or NULL");
        return result;
    }

    /* Initialise default parameters */
    s_fparams *params = init_def_fparams();
    if (!params) {
        result->error_code = 2;
        snprintf(result->error_message, sizeof(result->error_message),
                 "Failed to initialise fpocket parameters");
        return result;
    }

    strncpy(params->pdb_path, pdb_path, M_MAX_PDB_NAME_LEN - 1);
    params->pdb_path[M_MAX_PDB_NAME_LEN - 1] = '\0';
    params->fpocket_running = 1;
    params->db_run = 1;  /* suppress stdout noise */

    /* Apply user overrides (0 = use default) */
    if (min_alpha_size > 0) params->asph_min_size = min_alpha_size;
    if (max_alpha_size > 0) params->asph_max_size = max_alpha_size;
    if (min_spheres_per_pocket > 0) params->min_pock_nb_asph = min_spheres_per_pocket;
    if (clustering_distance > 0) params->clust_max_dist = clustering_distance;
    if (mc_iterations > 0) params->nb_mcv_iter = mc_iterations;

    /* Open and read PDB */
    s_pdb *pdb = open_and_read_pdb(pdb_path, 0, params);
    s_pdb *pdb_w_lig = open_and_read_pdb(pdb_path, 1, params);

    if (!pdb || !pdb_w_lig) {
        result->error_code = 3;
        snprintf(result->error_message, sizeof(result->error_message),
                 "Failed to read PDB file: %s", pdb_path);
        if (pdb) free_pdb_atoms(pdb);
        if (pdb_w_lig) free_pdb_atoms(pdb_w_lig);
        free_fparams(params);
        free_all();
        return result;
    }

    create_coord_grid(pdb);

    /* Run pocket detection */
    c_lst_pockets *pockets = search_pocket(pdb, params, pdb_w_lig);

    if (!pockets || pockets->n_pockets == 0) {
        result->num_pockets = 0;
        result->pockets = NULL;
        result->error_code = 0;
        if (pockets) c_lst_pocket_free(pockets);
        free_pdb_atoms(pdb);
        free_pdb_atoms(pdb_w_lig);
        free_fparams(params);
        free_all();
        return result;
    }

    /* Convert linked list to flat array */
    int n = (int)pockets->n_pockets;
    result->num_pockets = n;
    result->pockets = calloc(n, sizeof(fpocket_pocket_result));

    if (!result->pockets) {
        result->error_code = 4;
        snprintf(result->error_message, sizeof(result->error_message),
                 "Memory allocation failed for %d pockets", n);
        c_lst_pocket_free(pockets);
        free_pdb_atoms(pdb);
        free_pdb_atoms(pdb_w_lig);
        free_fparams(params);
        free_all();
        return result;
    }

    node_pocket *cur = pockets->first;
    for (int i = 0; i < n && cur; i++, cur = cur->next) {
        s_pocket *p = cur->pocket;
        fpocket_pocket_result *r = &result->pockets[i];

        r->rank = p->rank;
        r->score = p->score;
        r->num_alpha_spheres = p->nAlphaApol + p->nAlphaPol;
        r->num_apolar_alpha_spheres = p->nAlphaApol;
        r->num_polar_alpha_spheres = p->nAlphaPol;
        r->bary_x = p->bary[0];
        r->bary_y = p->bary[1];
        r->bary_z = p->bary[2];

        if (p->pdesc) {
            s_desc *d = p->pdesc;
            r->druggability_score = d->drug_score;
            r->volume = d->volume;
            r->convex_hull_volume = d->convex_hull_volume;
            r->hydrophobicity_score = d->hydrophobicity_score;
            r->polarity_score = (float)d->polarity_score;
            r->charge_score = (float)d->charge_score;
            r->prop_polar_atoms = d->prop_polar_atm;
            r->mean_alpha_sphere_radius = d->mean_asph_ray;
            r->alpha_sphere_density = d->as_density;
            r->mean_local_hydrophobic_density = d->mean_loc_hyd_dens;
            r->apolar_alpha_sphere_proportion = d->apolar_asphere_prop;
            r->max_alpha_sphere_distance = d->as_max_dst;
            r->flexibility = d->flex;
            r->surf_vdw = d->surf_vdw;
            r->surf_pol_vdw = d->surf_pol_vdw;
            r->surf_apol_vdw = d->surf_apol_vdw;
            r->inter_chain = d->interChain;
            r->n_abpa = d->n_abpa;
            strncpy(r->chain1, d->nameChain1, sizeof(r->chain1) - 1);
            strncpy(r->chain2, d->nameChain2, sizeof(r->chain2) - 1);
            strncpy(r->lig_tag, d->ligTag, sizeof(r->lig_tag) - 1);
            memcpy(r->aa_composition, d->aa_compo, sizeof(d->aa_compo));
        }
    }

    result->error_code = 0;
    c_lst_pocket_free(pockets);
    free_pdb_atoms(pdb);
    free_pdb_atoms(pdb_w_lig);
    free_fparams(params);
    free_all();

    return result;
}


void fpocket_free_result(fpocket_result *result) {
    if (!result) return;
    if (result->pockets) free(result->pockets);
    free(result);
}
