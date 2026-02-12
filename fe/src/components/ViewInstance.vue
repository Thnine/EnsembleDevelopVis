<template>
    <div>
        <GridVis1 v-if="drawData.type=='GridVis1'" :drawData="drawData"></GridVis1>
    </div>

</template>

<script>

import GridVis1 from '@/components/GridVis1.vue';

export default {
    name:"ViewInstance",
    props:['drawData'],
    components:{
        GridVis1
    },
    data(){
        return{

        }
    }

}
</script>

<style>

</style>