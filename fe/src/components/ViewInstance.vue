<template>
    <div>
        <GridVis1 v-if="drawData.type=='GridVis1'" :drawData="drawData"></GridVis1>
        <GridVis2 v-if="drawData.type=='GridVis2'" :drawData="drawData"></GridVis2>

    </div>

</template>

<script>

import GridVis1 from '@/components/GridVis1.vue';
import GridVis2 from '@/components/GridVis2.vue';

export default {
    name:"ViewInstance",
    props:['drawData'],
    components:{
        GridVis1,
        GridVis2
    },
    data(){
        return{

        }
    }

}
</script>

<style>

</style>